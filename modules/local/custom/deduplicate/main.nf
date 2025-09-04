// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

// FIXME: Eigen3 version is hard-coded because there is no CLI way to extract it

process Deduplicate {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "scottnortonphd/msdllcpapers/dragon:latest"

    input:
    tuple val(meta), path(annotatedBam), path(barcodes)
    path annofile
    val ambiguous

    output:
    tuple val(meta), path("${annotatedBam.baseName}_MarkDups.bam"), path(barcodes), emit: bam
    tuple val(meta), path("DRAGoN.out.${prefix}"), emit: counts
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    args = task.ext.args ?: ''
    """
deduplicate \\
    "${annotatedBam}" \\
    ${barcodes} \\
    ${annofile} \\
    "DRAGoN.out.${prefix}" \\
    --cpus ${task.cpus} \\
    --mmap-strategies ${ambiguous} \\
    ${args}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
    cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
    pkg-config: \$(pkg-config --version)
    zlib: \$(pkg-config --modversion zlib)
    bamtools: \$(pkg-config --modversion bamtools-1)
    eigen: 3.2
END_VERSIONS
"""

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    """
touch ${annotatedBam.baseName}_MarkDups.bam
mkdir DRAGoN.out.${prefix}
touch DRAGoN.out.${prefix}/{matrix.mtx,features.tsv,barcodes.tsv}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
    cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
    pkg-config: \$(pkg-config --version)
    zlib: \$(pkg-config --modversion zlib)
    bamtools: \$(pkg-config --modversion bamtools-1)
    eigen: 3.2
END_VERSIONS
"""
}
