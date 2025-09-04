// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process BuildStarIndex {
    label 'process_high'
    label 'process_high_memory'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'}"

    input:
    path genomefa
    path annofile
    val genome

    output:
    path "*", emit: index
    path 'versions.yml', emit: versions

    script:
    args = task.ext.args ?: ''
    annofile_arg = annofile ? "--sjdbGTFfile \"${annofile}\"" : ''
    true_fasta = file(genome.fasta).toUriString()
    if (annofile) {
        true_gtf = file(genome.gtf).toUriString()
        sed_command = "sed -i -r 's#${genomefa}#${true_fasta}#g; s#${annofile}#${true_gtf}#g' genomeParameters.txt"
    }
    else {
        sed_command = "sed -i -r 's#${genomefa}#${true_fasta}#g' genomeParameters.txt"
    }
    """
source <(guess-star-genome-params.awk "${genomefa}")
STAR --runMode genomeGenerate \\
    --runThreadN ${task.cpus} \\
    --genomeFastaFiles "${genomefa}" \\
    --genomeDir . \\
    ${annofile_arg} \\
    --genomeSAindexNbases \$genomeSAindexNbases \\
    --genomeChrBinNbits \$genomeChrBinNbits \\
    --limitGenomeGenerateRAM ${task.memory.bytes - 524288} \\
    ${args}
${sed_command}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
END_VERSIONS
"""

    stub:
    true_fasta = file(genome.fasta).toUriString()
    true_gtf = genome.gtf ? file(genome.gtf).toUriString() : '-'
    """
touch \\
    chrLength.txt \\
    chrNameLength.txt \\
    chrName.txt \\
    chrStart.txt \\
    exonGeTrInfo.tab \\
    exonInfo.tab \\
    geneInfo.tab \\
    Genome \\
    genomeParameters.txt \\
    Log.out \\
    SA \\
    SAindex \\
    sjdbInfo.txt \\
    sjdbList.fromGTF.out.tab \\
    sjdbList.out.tab \\
    transcriptInfo.tab
cat <<-END_PARAMS > genomeParameters.txt
versionGenome\t2.7.4a
genomeFastaFiles\t${true_fasta}
sjdbGTFfile\t${true_gtf}
END_PARAMS
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
END_VERSIONS
"""
}
