
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {CountUnmappedReads} from "../../../modules/local/custom/count_unmapped_reads/main.nf"
include {NormalizeCountsWF} from "../../../subworkflows/local/NormalizeCounts/main.nf"
include {STAR_ALIGN} from '../../../modules/local/star/align/main'
include {GenerateReport} from '../../../modules/local/custom/merge_solo/main'

workflow STARsolo {
    take:
        reads_in
        demux_log
        staridx
        annofile
        downsample_report
        val_ambiguous
        val_adjust_effective_lengths

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()

        STAR_ALIGN(
            reads_in,
            staridx,
            annofile,
            true
        )
        versions = versions.mix(STAR_ALIGN.out.versions.first())
        multiqc_files = multiqc_files.mix(STAR_ALIGN.out.log.map{ it -> it[1] })

        CountUnmappedReads (
            STAR_ALIGN.out.unmapped
                .map { meta, bam -> [meta.id, meta, bam] }
                .groupTuple()
                .map { _id, meta, bam ->
                    def x = [meta, bam]
                        .transpose()
                        .sort{ it -> it[0].bcidx }
                        .transpose()
                    x[0] = x[0][0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)}
                    x
                }
        )
        versions = versions.mix(CountUnmappedReads.out.versions.first())

        merged_output = STAR_ALIGN.out.solo
            .map { meta, outs, barcodes ->
                [meta.id, meta, outs, barcodes]
            }
            .groupTuple()
            .join(demux_log.map{ meta, json ->
                [ meta.id, json ]
            })
            .join(CountUnmappedReads.out.report.map{ meta, report ->
                [ meta.id, report ]
            })
            .join(downsample_report.map{ meta, report ->
                [ meta.id, report ]
            })
            .map { _id, meta, matrices, barcodes, json, unmap, downsample ->
                def x = [meta, matrices]
                    .transpose()
                    .sort{ it -> it[0].bcidx }
                    .transpose()
                x[0] = x[0][0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)}
                return x + [barcodes[0], json, unmap, downsample]
            }
        GenerateReport(
            merged_output,
            val_ambiguous
        )
        versions = versions.mix(GenerateReport.out.versions.first())

        counts_matrix = GenerateReport.out.report
            .map{ meta, _report, dense, _sparse ->
                [ meta, dense ]
            }
            .dump( tag: 'counts_matrix' )
        NormalizeCountsWF(
            STAR_ALIGN.out.bam,
            counts_matrix,
            annofile,
            val_adjust_effective_lengths
        )
        versions = versions.mix(NormalizeCountsWF.out.versions)

    emit:
        solo_bams = STAR_ALIGN.out.bam_sorted
        unmapped = STAR_ALIGN.out.unmapped
        counts_matrix
        normed_counts = NormalizeCountsWF.out.normed_counts
        versions
        multiqc_files
}
