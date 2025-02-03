#!/usr/bin/env nextflow

// Written by Scott Norton
// 12 Aug 2022

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {main_DRAGoN} from "./workflows/DRAGoN.nf"
include {main_STARsolo} from "./workflows/STARsolo.nf"

workflow {
    genome = params.private.GENOME_PATHS.get(params.IO.reference) ?: [:]
    genome.fasta = params.IO.fasta ?: genome.fasta
    genome.gtf = params.IO.gtf ?: genome.gtf
    genome.star = params.IO.star ?: genome.star
    genome.star_version = params.private.versionGenome
    ambiguous = (params.Dedup.ambiguous.split(",") + "Unique").toUnique().join(" ")

    if (params.pipeline == 'DRAGoN') {
        main_DRAGoN (
            params,
            genome,
            params.IO.fastqdir,
            params.IO.metadata,
            params.IO.bcldir,
            params.IO.sampleSheet,
            params.IO.sampleName,
            params.IO.barcodes,
            params.IO.contaminant,
            params.IO.no_prestage_index,
            params.Demux.dmuxn,
            params.PLATE_LAYOUT,
            params.Demux.data_quality_errors,
            params.IO.keep_bam,
            params.IO.outdir,
            params.Dedup.umi_mismatch,
            params.Demux.umilen,
            ambiguous,
            params.Count.adjust_effective_lengths,
        )
    } else if (params.pipeline == 'STARsolo') {
        main_STARsolo (
            params,
            genome,
            params.IO.fastqdir,
            params.IO.metadata,
            params.IO.bcldir,
            params.IO.sampleSheet,
            params.IO.sampleName,
            params.IO.barcodes,
            params.IO.contaminant,
            params.IO.no_prestage_index,
            params.Demux.dmuxn,
            params.PLATE_LAYOUT,
            params.Demux.data_quality_errors,
            params.IO.keep_bam,
            params.IO.outdir,
            ambiguous,
            params.Count.adjust_effective_lengths,
        )
    } else {
        def validPipelines = [
            'DRAGoN',
            'STARsolo'
        ]
        error "Invalid pipeline passed: --pipeline ${params.pipeline}\nChoose from: ${validPipelines.join(", ")}"
    }
}
