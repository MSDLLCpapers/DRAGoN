
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {Demultiplex} from '../../../modules/local/custom/demultiplex/main'
include {MergeDemuxJson} from '../../../modules/local/custom/merge_demux_json/main'
include {SAMTOOLS_SPLIT} from '../../../modules/local/samtools/split/main'
include {DownsampleBAM} from '../../../modules/local/custom/downsample/main'
include {Decontaminate} from "../../../subworkflows/local/Decontaminate/main.nf"
include {FastQC} from "../../../modules/local/fastqc/main.nf"

// Common workflow for preprocessing reads
// Demultiplexing, decontamination
workflow Preprocessing {
    take:
        reads_in
        val_dmuxn
        val_plate_layout
        val_data_quality_errors
        val_keep_bam
        val_contaminant

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()
        fastqc_in = Channel.empty()

        fastqc_in = fastqc_in.mix(reads_in.map{ meta, fastq, _barcodes ->
            [meta + [filetype: 'fastq', stage: 'untrimmed_fastqc'], fastq]
        })

        // Reinstate support for --Demux.dmuxn
        // LDDS-4848: Limit to 500M reads per job so as to keep file sizes low for transfer
        // Even though passing by: 0 will not split the file, I use this ternary to avoid the unnecessary overhead of reading and rewriting the file
        // FIXME: Race condition within nextflow causes truncated chunks to be dispatched
        if (val_dmuxn != 0 ) {
            reads_in
                .filter { _meta, fastq, _barcodes ->
                    fastq[0].scheme == 's3'
                }
                .first()
                .subscribe {
                    log.warn "BUG: Specifying nonzero --Demux.dmuxn with fastq from S3 can cause truncated input. Demultiplex will fail if this is encountered."
                }
            reads_in = reads_in
                .flatMap { meta, fastq, barcodes ->
                    [
                        [ meta, fastq[0], "R1", barcodes ],
                        [ meta, fastq[1], "R2", barcodes ]
                    ]
                }
                .splitFastq(by: val_dmuxn, file: true, elem: 1)
                .map{ meta, fastq, end, barcodes ->
                    [ meta.id2, fastq.baseName.tokenize('_')[-1], meta, fastq, end, barcodes ]
                }
                .groupTuple(by: [0, 1])
                .map{ _id2, chunk, meta_l, fastq_l, end_l, barcodes_l ->
                    assert meta_l.size() == 2
                    def fastq = [ fastq_l, end_l ]
                        .transpose()
                        .sort{ it -> it[1] }
                        .transpose()
                        .getAt(0)
                    return [ meta_l[0] + [ id2 : meta_l[0].id2 + "_${chunk}" ], fastq, barcodes_l[0] ]
                }
            reads_in.dump(tag: 'after_split')
        }
        Demultiplex(reads_in)
        versions = versions.mix(Demultiplex.out.versions.first())

        Decontaminate (
            Demultiplex.out.bam,
            val_contaminant
        )
        versions = versions.mix(Decontaminate.out.versions)
        multiqc_files = multiqc_files.mix(Decontaminate.out.multiqc_files)
        contaminated = Channel.empty()
            .mix(
                Decontaminate.out.contaminated,
                Demultiplex.out.qcfail
            )

        fastqc_in = fastqc_in.mix(Decontaminate.out.decontaminated.map{ meta, bam, _barcodes ->
            [meta + [filetype: 'bam', stage: 'trimmed_fastqc'], bam]
        })
        fastqc_in.dump (tag: 'fastqc_in')
        FastQC (fastqc_in, false)
        versions = versions.mix(FastQC.out.versions.first())
        multiqc_files = multiqc_files.mix(FastQC.out.zip)

        SAMTOOLS_SPLIT(Decontaminate.out.decontaminated)
        versions = versions.mix(SAMTOOLS_SPLIT.out.versions.first())

        // DSSI-111: QC checkpoint
        // DSSI-451: Add key "plate_map"
        MergeDemuxJson(
            Demultiplex.out.log
                .map{meta, json, barcodes ->
                    [meta.id, meta, json, barcodes]
                }
                .groupTuple()
                .map{_id, meta, json, barcodes ->
                    def x = [meta, json]
                        .transpose()
                        .sort{ it -> it[0].id2 }
                        .transpose()
                    x[0] = x[0][0].findAll{k, _v -> k != 'id2'}
                    return x + [barcodes[0]]
                },
            val_plate_layout,
            val_data_quality_errors
        )
        versions = versions.mix(MergeDemuxJson.out.versions.first())

        // This ugly mess is to create a hard checkpoint for error detection,
        // and arrest the pipeline if QC issues are to be treated as fatal.
        // It also parses the result of SAMTOOLS_SPLIT to separate the channel
        // by well.

        SAMTOOLS_SPLIT.out.bam
            .map{meta, bam, barcodes ->
                [meta.id, meta, bam, barcodes]
            }
            .groupTuple()
            .combine(MergeDemuxJson.out.merged.map{meta, json, stdout ->
                [meta.id, json, stdout]
            }, by: 0)
            .tap { ch_split_json_stdout }
            .flatMap{ _id, meta, bam, barcodes, json, stdout ->
                if (stdout) {
                    if (val_data_quality_errors != "ignore") {
                        log.warn "${meta[0].id}: ${stdout}"
                    }
                }
                bam
                    .flatten()
                    .collect{fl ->
                        return (fl.baseName =~ /^(\S+)_unmap_(\w{4})$/)
                            .collect{ m ->
                                def bcidx = m[2]
                                def bcname
                                if (bcidx != 'noBC') {
                                    bcname = barcodes[0]
                                        .readLines()
                                        .getAt(bcidx as Integer)
                                        .strip()
                                        .tokenize("\t")
                                        .getAt(1)
                                } else {
                                    bcname = 'noBC'
                                }
                                [
                                    meta[0] + [
                                        id2 : m[1],
                                        'bcidx': bcidx,
                                        'bcname': bcname,
                                    ],
                                    fl,
                                    barcodes[0],
                                    json
                                ]
                            }
                    }
                        .flatten()
                        .collate(4)
            }
            .dump(tag: 'after_flatMap')
            .filter{meta, _bam, _barcodes, _demux_json ->
                val_keep_bam || meta.bcidx != "noBC"
            }
            .tap { downsample_in }
            .map { it -> [false] + it }
            .set { ch_split_bam }

        ch_split_json_stdout
            .filter { _id, _meta, _bam, _barcodes, _json, stdout ->
                stdout
            }
            .map { _id, meta, _bam, _barcodes, _json, stdout ->
                if (val_data_quality_errors != "ignore") {
                    log.warn "${meta[0].id}: ${stdout}"
                }
                return true
            }
            .ifEmpty { false }
            .first()
            .set { ch_any_errors }

        ch_split_bam
            .combine(ch_any_errors, by: 0)
            .map { _x, meta, bam, barcodes, demux_json ->
                [ meta, bam, barcodes, demux_json ]
            }
            .tap { ch_downsample_in }
            .ifEmpty {
                if (val_data_quality_errors == "raise") {
                    log.error "Data quality problems detected. Workflow will now exit. See ${launchDir}/.nextflow.log for details."
                } else if (val_data_quality_errors == "warn") {
                    log.warn "Data quality problems detected. Workflow will proceed. See ${launchDir}/.nextflow.log for details."
                }
            }

        if (val_data_quality_errors == 'raise') {
            downsample_in = ch_downsample_in
        }

        downsample_in
            .branch { meta, _bam, _barcodes, _demux ->
                nobc: meta.bcidx == 'noBC'
                has_bc: true
            }
            .set { ch_downsample_branch }

        // Run this constitutively to generate the report and remove qcfail reads
        DownsampleBAM(ch_downsample_branch.has_bc)
        versions = versions.mix(DownsampleBAM.out.versions.first())

        DownsampleBAM.out.bam
            .map{ meta, bam, barcodes ->
                [ meta.id, meta.bcidx, meta, bam, barcodes ]
            }
            .groupTuple(by: [0, 1])
            .map{ _id, _bcidx, meta, bam, barcodes ->
                def x = [meta, bam]
                    .transpose()
                    .sort{ it -> it[0].id2 }
                    .transpose()
                x[0] = x[0][0].findAll{k, _v -> k != 'id2'}
                return x + [ barcodes[0] ]
            }
            .set {split_bam}

        DownsampleBAM.out.report
            .collectFile (keepHeader: true, skip: 1, sort: { content -> content.tokenize('\n').last().tokenize('\t')[2] as Integer }) { meta, fl ->
                [ "${meta.id2}.downsample.txt", fl.text ]
            }
            .map { it ->
                [ (it.name =~ /^(\S+)\.downsample\.txt/)[0][1], it ]
            }
            .join(DownsampleBAM.out.report.map {meta, _fl -> [meta.id2, meta]})
            .map { _id2, fl, meta ->
                [ meta.id, meta, fl ]
            }
            .groupTuple()
            .map { _id, meta, fl ->
                def x = [meta, fl]
                    .transpose()
                    .sort{ it -> it[0].id2 }
                    .transpose()
                x[0] = x[0][0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)}
                x
            }
            .set {downsample_report}

        ch_downsample_branch
            .nobc
            .map{meta, bam, _barcodes, _demux ->
                [meta, bam]
            }
            .mix(DownsampleBAM.out.discard)
            .set {discarded_bam}

    // FIXME: workflow is null inside onComplete handler
    def myWorkflow
    myWorkflow = workflow
    workflow.onComplete {
        ch_any_errors.subscribe { flag ->
            if (flag) {
                if (val_data_quality_errors == "raise") {
                    error "Reporting failed execution due to data quality issues."
                } else if (val_data_quality_errors == "warn") {
                    log.warn "Workflow completed. Some wells had data quality issues. See ${myWorkflow.launchDir}/.nextflow.log for details."
                }
            }
        }
    }

    emit:
        split_bam
        decontam = Decontaminate.out.decontaminated
        contaminated
        demux_log = MergeDemuxJson.out.merged.map{ meta, json, _stdout ->
            [ meta, json ]
        }
        downsample_report
        discarded_bam
        versions
        multiqc_files
        any_errors = ch_any_errors
}

