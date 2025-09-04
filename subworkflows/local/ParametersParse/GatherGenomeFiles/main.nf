// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include { BCL2FASTQ      } from '../../../../modules/local/bcl2fastq/main'
include { BuildStarIndex } from '../../../../modules/local/star/genomegenerate/main'
include { PIGZ as Unpigz } from '../../../../modules/local/pigz/main'
include { StageStarIndex } from '../../../../modules/local/custom/stage_star_index/main'

workflow GatherGenomeFiles {
    take:
    val_genome
    val_no_prestage_index

    main:
    versions = Channel.empty()
    multiqc_files = Channel.empty()
    pigz_files_in = Channel.empty()

    // STAR reference always needs to be defined
    // Deprecated: Only require STAR index if plan to align with STAR
    staridx = val_genome.star
    if (staridx == null) {
        error("Missing required argument: --IO.star <DIR> (hint: you can provide an alias (GRCh38, mm10) via --IO.reference)")
    }

    def need_build_star_index
    if (file("${staridx}/genomeParameters.txt").exists()) {
        need_build_star_index = false
        // To be deprecated: automatic inference of fasta and gtf files from STAR index
        // Will use the genome map or --IO.fasta, --IO.gtf, see below
        (genomefa, annofile) = ValidateStarIndex(val_genome)
        staridx = Channel.of(file(staridx, type: 'dir').listFiles()).collect()
        if (file(workDir).scheme != 's3' && !val_no_prestage_index) {
            StageStarIndex(staridx)
            staridx = StageStarIndex.out.files
            versions = versions.mix(StageStarIndex.out.versions)
        }
    }
    else {
        need_build_star_index = true
        genomefa = Channel.value(file(val_genome.fasta, checkIfExists: true))
        annofile = Channel.value(file(val_genome.gtf, checkIfExists: true))
    }

    // Fasta and GTF files may be gzipped
    pigz_files_in = pigz_files_in.mix(
        genomefa.filter { fn -> fn.extension == 'gz' }.map { it -> ['genomefa', it] },
        annofile.filter { fn -> fn.extension == 'gz' }.map { it -> ['annofile', it] },
    )
    Unpigz(pigz_files_in)
    versions = versions.mix(Unpigz.out.versions.first())
    pigz_files_out = Unpigz.out.result
        .view { meta, fp ->
            "Uncompressed ${meta} output at: ${fp.toUriString()}"
        }
        .branch { meta, _fp ->
            genomefa: meta == 'genomefa'
            annofile: meta == 'annofile'
            unknown: true
        }
    genomefa = genomefa
        .filter { fn -> fn.extension != 'gz' }
        .mix(pigz_files_out.genomefa.map { _meta, fp -> fp })
        .collect()
    annofile = annofile
        .filter { fn -> fn.extension != 'gz' }
        .mix(pigz_files_out.annofile.map { _meta, fp -> fp })
        .collect()

    // Call STAR genomeGenerate as needed
    if (need_build_star_index) {
        log.info("STAR index at ${staridx} not found, will build index there")
        BuildStarIndex(
            genomefa,
            annofile,
            val_genome,
        )
        staridx = BuildStarIndex.out.index
        versions = versions.mix(BuildStarIndex.out.versions)
    }

    emit:
    staridx
    genomefa
    annofile
    versions
    multiqc_files
}

// =============================
// STAR index
// =============================

// Make sure the STAR index has the correct version, and grab the GTF file from it
def ValidateStarIndex(genome) {
    log.warn("DeprecationWarning: automatic inference of the genome fasta and annotation paths from the STAR index path is deprecated and will be removed in a future release")
    def genomeParameters = file("${genome.star}/genomeParameters.txt")
        .readLines()
        .findAll { line -> line != ~/^#/ }
        .collectEntries { line -> line.split('\t', 2) }

    if (genomeParameters.versionGenome != genome.star_version) {
        error("Expected version ${genome.star_version}, got ${genomeParameters.versionGenome}")
    }

    def annofile = []
    if (genomeParameters.getOrDefault('sjdbGTFfile', '-') != '-') {
        annofile = genomeParameters.sjdbGTFfile.split(' ')
    }
    annofile = annofile.findAll { fl -> file(fl).exists() }
    if (annofile.isEmpty()) {
        if (!genome.gtf) {
            error("Unable to deduce GTF from STAR index, and --IO.gtf not provided")
        }
        if (!file(genome.gtf).exists()) {
            error("${genome.gtf}: no such file or directory")
        }
        annofile = [genome.gtf]
    }

    def genomefa = genomeParameters.getOrDefault('genomeFastaFiles', '').split(' ').findAll { it -> it != '' && file(it).exists() }
    if (genomefa.isEmpty()) {
        if (!genome.fasta) {
            error("Unable to deduce reference fasta from STAR index, and --IO.fasta not provided")
        }
        if (!file(genome.fasta).exists()) {
            error("${genome.fasta}: no such file or directory")
        }
        genomefa = [genome.fasta]
    }

    return [
        Channel.value(genomefa.collect { fn -> file(fn, checkIfExists: true) }),
        Channel.value(annofile.collect { fn -> file(fn, checkIfExists: true) }),
    ]
}
