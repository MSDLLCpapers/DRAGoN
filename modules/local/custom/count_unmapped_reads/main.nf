#!/usr/bin/env nextflow

// Written by Scott Norton
// 07 May 2024

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process CountUnmappedReads {
    label 'process_single'
    input:
        tuple val(meta), path(unmappedBam)
    output:
        tuple val(meta), path("${meta.id}_unmapped.txt"), emit: report
        path 'versions.yml', emit: versions
    script:
"""
samtools cat -@ ${task.cpus} ${unmappedBam} | samtools view -F QCFAIL | grep -oE 'BI:Z:\\w+' | sort | uniq -c | sed -r 's/^\\s+([0-9]+)\\s+BI:Z:(\\w+)/\\2\\t\\1/g' > ${meta.id}_unmapped.txt
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
    stub:
"""
touch "${meta.id}_unmapped.txt"
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
"""
}
