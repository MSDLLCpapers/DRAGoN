
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include { BCL2FASTQ                        } from '../../../../modules/local/bcl2fastq/main'
include { MaybeRestoreFromS3               } from '../../../../modules/local/custom/s3_archive_restore/main'
include { ParseMetadataXLSX                } from '../../../../modules/local/custom/parse_metadata/main'

// =============================
// Fastq files
// =============================

def makeMeta(fastqFile) {
    def matches = fastqFile.name =~ /^((.+)_S\d+(?:_L\d+)?)_(R[12])(_.+)?\.(fastq|fq)(\.gz)?/
    matches.size() ? [ id: matches[0][2], id2: matches[0][1] + matches[0][4] ] : 0
}

workflow GatherInputFastq {
    take:
        val_fastqdir
        val_metadata
        val_bcldir
        val_sampleSheet
        val_sampleName
        val_barcodes

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()

        // Validate metadata and samplesheet
        sampleSheet_in = val_sampleSheet ? file(val_sampleSheet, checkIfExists: true) : []
        if (val_metadata != null) {
            metadata = Channel.fromPath(val_metadata)
            ParseMetadataXLSX(metadata, sampleSheet_in)
            sampleSheet = ParseMetadataXLSX.out.sampleSheet.first()
            versions = versions.mix(ParseMetadataXLSX.out.versions.first())
            samples_in = ParseMetadataXLSX.out.barcodes
                .flatten()
                .map{ fname ->
                    [ fname.simpleName.replace("_barcodes", ""), fname ]
                }
        } else {
            samples_in = Channel.of([val_sampleName, file(val_barcodes, checkIfExists: true)])
            sampleSheet = Channel.of(sampleSheet_in)
        }
        samples = samples_in.map{ sampleName, barcodes ->
            def lines = barcodes.readLines()
            def meta = [
                id: sampleName,
                numbc: lines.size,
                bclen: lines[0].tokenize("\t")[0].length()
            ]
            return [ meta.id, meta, barcodes ]
        }
            .dump(tag: 'samples')

        // Check if fastq need to be generated or retrieved from S# archive
        // Assume fastq come from BCL2FASTQ
        filesGlob = "${val_fastqdir}/**_R[12]_*.{fastq,fq}{,.gz}"
        if (file(filesGlob).isEmpty()) {
            try{
                bcldir = file(val_bcldir, type: "dir", checkIfExists: true)
            } catch (_all) {
                error "Fastq directory ${val_fastqdir} is empty, --IO.bcldir is required"
            }
            sampleSheet.filter{it -> it}.ifEmpty{
                error "Fastq directory ${val_fastqdir} is empty, --IO.sampleSheet is required"
            }
            log.info "Fastq directory ${val_fastqdir} is empty, will attempt to populate using bcl2fastq"
            BCL2FASTQ(bcldir, sampleSheet)
            fastq_files = BCL2FASTQ.out.fastq
                .flatMap{ it ->
                    it.collect { fastq ->
                        def meta = makeMeta(fastq)
                        [ meta, fastq ]
                    }
                }
                .groupTuple(sort: { it ->
                    it.name
                })
            versions = versions.mix(BCL2FASTQ.out.versions)
            multiqc_files = multiqc_files.mix(BCL2FASTQ.out.stats)
        } else {
            fastq_files_check = Channel.fromFilePairs(filesGlob) { file ->
                makeMeta(file)
            }
                .filter{meta, _fastq -> meta.id != 'Undetermined'}
                .branch{_meta, fastq ->
                    s3: fastq.findAll{ it -> it.scheme == 's3' }
                    rest: true
                }
            fastq_files_check.rest.dump(tag: 'maybe_glacial_rest')
            maybe_glacial = fastq_files_check.s3
                .map{meta, fastq -> [meta, fastq.collect{ it -> it.toUriString() }]}
                .ifEmpty{[null, []]}
                .toList()
                .map{collected ->
                    collected.transpose()
                }
                .dump(tag: 'maybe_glacial_s3')
            MaybeRestoreFromS3(maybe_glacial)
            fastq_files = MaybeRestoreFromS3.out.fastq
                .transpose()
                .mix(fastq_files_check.rest)
                .dump(tag: 'after_restore')
            versions = versions.mix(MaybeRestoreFromS3.out.versions)
        }

        // Match fastq files with metadata
        reads_in_tocross = fastq_files
            .filter{meta, _fileNames -> meta}
            .map{meta, fileNames ->
                [meta.id, meta, fileNames]  // put sampleName first for cross
            }
            .groupTuple()

        reads_in = samples
            .join(reads_in_tocross, remainder: true) // Map fastq files to barcodes
            .transpose()
            .filter{ it -> it.size() == 5 }
            .map{id, meta, barcodes, meta2, fastq -> // Rearrange to [meta, [R1, R2], barcodes]
                if (!fastq) {
                    error "unable to find fastq corresponding to sample: ${id}"
                }
                [ meta + meta2, fastq, barcodes ]
            }
            .filter{meta, _fastq, _barcodes ->
                val_sampleName ? val_sampleName.split(",").contains(meta.id) : true
            }
            .ifEmpty{
                error (val_sampleName ? "--IO.sampleName specified but none matches the fastq files" : "Failure matching metadata with fastq")
            }
            .dump (tag: 'reads_in')

    emit:
        reads_in // tuple val(meta), path(barcodes), path(reads)
        versions
        multiqc_files
}
