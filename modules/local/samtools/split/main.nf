
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process SAMTOOLS_SPLIT {
    label 'process_medium'

    input:
        tuple val(meta), path(mergedbam), path(barcodes)

    output:
        tuple val(meta), path("${prefix}*.bam"), path(barcodes), emit: bam
        path 'versions.yml', emit: versions

    script:
        prefix = task.ext.prefix ?: meta.id
        args = task.ext.args ?: ''
"""
samtools \\
    split \\
    -@ ${task.cpus} \\
    -f "${prefix}_%!.%." \\
    ${args} \\
    ${mergedbam}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: meta.id
"""
for i in {0..${meta.numbc - 1}}; do
    touch ${prefix}_\$(printf "%04d" \$i).bam
done
touch ${prefix}_noBC.bam
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
}
