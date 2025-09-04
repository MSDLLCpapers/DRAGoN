#!/usr/bin/env nextflow

// Written by Scott Norton
// 12 Aug 2022

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include { main_DRAGoN   } from "./workflows/DRAGoN.nf"
include { main_STARsolo } from "./workflows/STARsolo.nf"

workflow {
    genome = params.genomes.get(params.IO.reference) ?: [:]
    genome.fasta = params.IO.fasta ?: genome.fasta
    genome.gtf = params.IO.gtf ?: genome.gtf
    genome.star = params.IO.star ?: genome.star
    genome.star_version = params.Align.versionGenome
    ambiguous = (params.Dedup.ambiguous.split(",") + "Unique").toUnique().join(" ")

    if (params.pipeline == 'DRAGoN') {
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
    else if (params.pipeline == 'STARsolo') {
        main_STARsolo(
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
    }
    else {
        def validPipelines = [
            'DRAGoN',
            'STARsolo',
        ]
        error("Invalid pipeline passed: --pipeline ${params.pipeline}\nChoose from: ${validPipelines.join(", ")}")
    }

    def myParams
    def myWorkflow

    myParams = params
    myWorkflow = workflow

    workflow.onComplete {
        def emails = myParams.notifyEmails
        if (emails instanceof List) {
            emails = emails.join(",")
        }
        if (emails != '') {
            log.info("Sending out emails to notify: ${emails}")
            def msg
            if (myWorkflow.success) {
                msg = """\
                    ${myParams.pipeline} ${myWorkflow.sessionId} completed.
                    ---------------------------
                    Parameters:
                    ---------------------------
                    IO.metadata: ${myParams.IO.metadata}
                    IO.fastqdir: ${myParams.IO.fastqdir}
                    IO.outdir: ${myParams.IO.outdir}
                    ---------------------------
                    Completed at: ${myWorkflow.complete}
                    Duration    : ${myWorkflow.duration}
                    Success     : ${myWorkflow.success}
                    workDir     : ${myWorkflow.workDir}
                    exit status : ${myWorkflow.exitStatus}
                    """.stripIndent()
            }
            else {
                msg = """\
                    ${myParams.pipeline} ${myWorkflow.sessionId} failed:
                    ---------------------------
                    Parameters:
                    ---------------------------
                    IO.metadata: ${myParams.IO.metadata}
                    IO.fastqdir: ${myParams.IO.fastqdir}
                    IO.outdir: ${myParams.IO.outdir}
                    Launch directory: ${myWorkflow.launchDir}
                    The launch directory contains full log file (.nextflow.log) and the HTML report of the run, which could be useful for troubleshooting.
                    ---------------------------
                    Completed at: ${myWorkflow.complete}
                    Duration    : ${myWorkflow.duration}
                    Success     : ${myWorkflow.success}
                    workDir     : ${myWorkflow.workDir}
                    exit status : ${myWorkflow.exitStatus}
                    ---------------------------
                    Error report:
                    ---------------------------
                    ${myWorkflow.errorReport}
                    """.stripIndent()
            }
            sendMail(
                to: emails,
                subject: "DRUG-seq (using ${myParams.pipeline}): ${myWorkflow.sessionId}",
                body: msg,
            )
        }
    }
}
