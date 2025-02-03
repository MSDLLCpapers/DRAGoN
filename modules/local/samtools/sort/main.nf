
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process SAMTOOLS_SORT {
    label 'process_medium'
    input:
        tuple val(meta), path(bam), path(barcodes)
    output:
        tuple val(meta), path("${prefix}.bam"), path(barcodes), emit: bam
        path 'versions.yml', emit: versions
    script:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
        args = task.ext.args ?: ''
        if (bam.baseName == prefix) {
            error "output file will have the same name as the input"
        }
"""
samtools \\
    sort \\
    -@ ${task.cpus} \\
    $args \\
    -o ${prefix}.bam \\
    $bam
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
"""
touch ${prefix}.bam
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
}
