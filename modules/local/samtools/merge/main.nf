// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process SAMTOOLS_MERGE {
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'
        : 'biocontainers/samtools:1.20--h50ea8bc_0'}"

    input:
    tuple val(meta), path(inbam, arity: '1..*')

    output:
    path "${prefix}.bam", emit: bam
    path "${prefix}.bam.bai", emit: bai
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (inbam.contains("${prefix}.bam")) {
        error("input and output filenames should not be the same")
    }
    """
samtools merge -f -c -p -@ ${task.cpus} -O BAM ${prefix}.bam ${inbam}
samtools index -@ ${task.cpus} ${prefix}.bam
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
touch ${prefix}.bam ${prefix}.bam.bai flagstat.json
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
}
