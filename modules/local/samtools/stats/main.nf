
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process SAMTOOLS_STATS {
    label 'process_medium'
    input:
        tuple val(meta), path(inbam), path(barcodes)
    output:
        tuple val(meta), path("${prefix}.stats"), emit: stats
        path 'versions.yml', emit: versions
    script:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
        args = task.ext.args ?: ''
"""
samtools \\
    stats \\
    -@ ${task.cpus} \\
    $args \\
    $inbam \\
    > ${prefix}.stats
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
    stub:
"""
export FRAGMENT_LENGTH=1048576
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
}
