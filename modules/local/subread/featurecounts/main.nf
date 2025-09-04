// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process FeatureCounts {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_2'
        : 'biocontainers/subread:2.0.6--he4a0461_2'}"

    input:
    tuple val(meta), path(alignedBam), path(barcodes)
    path annofile

    output:
    tuple val(meta), path("*featureCounts.bam"), path(barcodes), emit: bam
    tuple val(meta), path("${prefix}.out.summary"), emit: stats
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    args = task.ext.args ?: ''
    """
featureCounts \\
    ${args} \\
    -a "${annofile}" -F GTF \\
    -M -O --fraction -T 1 \\
    -R BAM -o "${prefix}.out" \\
    "${alignedBam}"
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    subread: \$(featureCounts -v 2>&1 | grep -v '^\$' | sed 's/.* v//g')
END_VERSIONS
"""

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    """
touch ${alignedBam}.featureCounts.bam ${prefix}.out.summary
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    subread: \$(featureCounts -v 2>&1 | grep -v '^\$' | sed 's/.* v//g')
END_VERSIONS
"""
}
