
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include { getWorkflowVersion } from '../../../subworkflows/local/util_drugseq/main'
include { GatherInputFastq } from '../../../subworkflows/local/ParametersParse/GatherInputFastq/main'
include { GatherGenomeFiles } from '../../../subworkflows/local/ParametersParse/GatherGenomeFiles/main'
include { validateParameters
          paramsSummaryMap
          paramsSummaryLog             } from 'plugin/nf-schema'

//
// Inputs validation workflow
//

workflow ParametersParse {
    take:
        val_params_copy
        val_genome
        val_fastqdir
        val_metadata
        val_bcldir
        val_sampleSheet
        val_sampleName
        val_barcodes
        val_no_prestage_index

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()

        log.info"""
        =================================================
        DRUG-seq Pipeline ${getWorkflowVersion()}
        =================================================
        """.stripIndent()

        try {
            file("${projectDir}/bin/.commit").withWriter{ wr -> wr << getWorkflowVersion()}
        } catch (_all) {}

        validateParameters ()

        log.info paramsSummaryLog (workflow)

        GatherInputFastq (
            val_fastqdir,
            val_metadata,
            val_bcldir,
            val_sampleSheet,
            val_sampleName,
            val_barcodes,
        )
        versions = versions.mix(GatherInputFastq.out.versions)
        multiqc_files = multiqc_files.mix(GatherInputFastq.out.multiqc_files)

        GatherGenomeFiles (
            val_genome,
            val_no_prestage_index
        )
        versions = versions.mix(GatherGenomeFiles.out.versions)
        multiqc_files = multiqc_files.mix(GatherGenomeFiles.out.multiqc_files)

        dumpParams(
            val_params_copy,
            GatherGenomeFiles.out.genomefa,
            GatherGenomeFiles.out.annofile,
            val_genome
        )

    emit:
        reads_in = GatherInputFastq.out.reads_in
        staridx = GatherGenomeFiles.out.staridx
        genomefa = GatherGenomeFiles.out.genomefa
        annofile = GatherGenomeFiles.out.annofile
        versions
        multiqc_files
}

// =============================
// Main entry
// =============================

// Dump parameters to yaml. Useful for resubmitting a failed run
def dumpParams(params, genomefa, annofile, genome) {
    Channel.of(['star', file(genome.star)])
        .mix(
            genomefa.map{it -> ['fasta', it]},
            annofile.map{it -> ['gtf', it]}
        )
        .groupTuple()
        .map{key, value ->
            value = value.flatten().collect{fname -> fname.toUriString()}
            return [ key, (value.size() == 1 ? value[0] : value) ]
        }
        .toList()
        .collectFile (name: "params.run.yaml", storeDir: params.IO.outdir) { tup ->
            def to_dump = [:] << params
            to_dump.remove("private")
            to_dump.remove("aws")
            to_dump['pipelineVersion'] = getWorkflowVersion()
            to_dump['IO']['outdir'] = to_dump['IO']['outdir'] ? file(to_dump['IO']['outdir']).toUriString() : null
            to_dump['IO']['metadata'] = to_dump['IO']['metadata'] ? file(to_dump['IO']['metadata']).toUriString() : null
            to_dump['IO']['fastqdir'] = to_dump['IO']['fastqdir'] ? file(to_dump['IO']['fastqdir']).toUriString() : null
            to_dump['IO']['bcldir'] = to_dump['IO']['bcldir'] ? file(to_dump['IO']['bcldir']).toUriString() : null
            to_dump['IO']['sampleSheet'] = to_dump['IO']['sampleSheet'] ? file(to_dump['IO']['sampleSheet']).toUriString() : null
            to_dump['IO']['contaminant'] = to_dump['IO']['contaminant'] ? file(to_dump['IO']['contaminant']).toUriString() : null
            to_dump['IO'] += tup.collectEntries()
            def options = new org.yaml.snakeyaml.DumperOptions()
            options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
            def yaml = new org.yaml.snakeyaml.Yaml(options)
            yaml.dump(to_dump)
        }
    if (params.IO.metadata) {
        file(params.IO.metadata).copyTo("${params.IO.outdir}/${file(params.IO.metadata).name}")
    }
    if (params.IO.sampleSheet) {
        file(params.IO.sampleSheet).copyTo("${params.IO.outdir}/${file(params.IO.sampleSheet).name}")
    }
}
