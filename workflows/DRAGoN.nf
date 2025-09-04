#!/usr/bin/env nextflow
include { softwareVersionsToYAML } from '../subworkflows/local/util_drugseq/main'
include { ParametersParse        } from "../subworkflows/local/ParametersParse/main.nf"
include { Preprocessing          } from "../subworkflows/local/Preprocessing/main.nf"
include { DRAGoN                 } from "../subworkflows/local/DRAGoN/main.nf"
include { Postprocessing         } from "../subworkflows/local/Postprocessing/main.nf"
include { MultiQC                } from "../modules/local/multiqc/main"


// Written by Scott Norton
// 07 Oct 2022

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

params.pipeline = "DRAGoN"

workflow {
    genome = params.genomes.get(params.IO.reference) ?: [:]
    genome.fasta = params.IO.fasta ?: genome.fasta
    genome.gtf = params.IO.gtf ?: genome.gtf
    genome.star = params.IO.star ?: genome.star
    genome.star_version = params.Align.versionGenome
    ambiguous = (params.Dedup.ambiguous.split(",") + "Unique").toUnique().join(" ")

    main_DRAGoN(
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
}

workflow main_DRAGoN {
    take:
    val_params_copy
    val_genome
    val_fastqdir
    val_metadata
    val_bcldir
    val_sampleSheet
    val_sampleName
    val_barcodes
    val_contaminant
    val_no_prestage_index
    val_dmuxn
    val_plate_layout
    val_data_quality_errors
    val_keep_bam
    val_outdir
    val_umi_mismatch
    val_umilen
    val_ambiguous
    val_adjust_effective_lengths

    main:
    // Initialize
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Resolve metadata and genome references
    ParametersParse(
        val_params_copy,
        val_genome,
        val_fastqdir,
        val_metadata,
        val_bcldir,
        val_sampleSheet,
        val_sampleName,
        val_barcodes,
        val_no_prestage_index,
    )
    ch_versions = ch_versions.mix(ParametersParse.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ParametersParse.out.multiqc_files)

    // Filtering, demultiplexing, and QC
    if (val_umi_mismatch * 3 >= val_umilen) {
        error("Allowing too many mismatches in the UMI sequence will impede PCR duplicate detection. Cowardly refusing to proceed.")
    }
    Preprocessing(
        ParametersParse.out.reads_in,
        val_dmuxn,
        val_plate_layout,
        val_data_quality_errors,
        val_keep_bam,
        val_contaminant,
    )
    ch_versions = ch_versions.mix(Preprocessing.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(Preprocessing.out.multiqc_files)

    // Alignment, deduplication, and counting
    DRAGoN(
        Preprocessing.out.split_bam,
        Preprocessing.out.demux_log,
        ParametersParse.out.staridx,
        ParametersParse.out.annofile,
        Preprocessing.out.downsample_report,
        val_ambiguous,
        val_adjust_effective_lengths,
        val_keep_bam,
    )
    ch_versions = ch_versions.mix(DRAGoN.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DRAGoN.out.multiqc_files)

    // Generate final BAM file
    Postprocessing(
        DRAGoN.out.marked_bam,
        Preprocessing.out.contaminated,
        DRAGoN.out.unmapped,
        Preprocessing.out.discarded_bam,
        val_keep_bam,
    )
    ch_versions = ch_versions.mix(Postprocessing.out.versions)

    // Dump version information
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${val_outdir}",
            name: "mrk_dependent_software_mqc_versions.yml",
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    // Final multiqc report
    MultiQC(
        ch_multiqc_files.collect(),
        "${projectDir}/assets/multiqc_config.yaml",
    )
}
