// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

// FIXME: Eigen3 version is hard-coded because there is no CLI way to extract it
// FIXME: Boost version is hard-coded because there is no CLI way to extract it

// Demultiplex the reads. STARsolo has its own
// demultiplexer but it's limited to 0 or 1 mismatches.
// cutadapt also has its own demultiplexer but it does not seem to work properly
// This approach accepts 2 or more mismatches and also
// splits the input into separate files for input.
process Demultiplex {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "scottnortonphd/msdllcpapers/dragon:latest"

    input:
    tuple val(meta), path(fastq), path(barcodes)

    output:
    tuple val(meta), path("${prefix}.bam"), path(barcodes), emit: bam
    tuple val(meta), path("*.demux.json"), path(barcodes), emit: log
    tuple val(meta), path("${prefix}_qcFail.bam"), optional: true, emit: qcfail
    path 'versions.yml', emit: versions

    script:
    prefix = "${task.ext.prefix ?: meta.id}"
    args = task.ext.args ?: ''
    """
demultiplex \\
    ${fastq} \\
    "${barcodes}" \\
    "${prefix}" \\
    --cpus ${task.cpus} \\
    ${args}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
    cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
    pkg-config: \$(pkg-config --version)
    zlib: \$(pkg-config --modversion zlib)
    nlohmann_json: \$(pkg-config --modversion nlohmann_json)
    kseq++: \$(pkg-config --modversion kseq++)
    bamtools: \$(pkg-config --modversion bamtools-1)
    eigen: 3.2
    boost: 1.84.0
END_VERSIONS
"""

    stub:
    prefix = "${task.ext.prefix ?: meta.id}"
    """
touch ${prefix}.bam ${prefix}_qcFail.bam ${prefix}.demux.json
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
    cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
    pkg-config: \$(pkg-config --version)
    zlib: \$(pkg-config --modversion zlib)
    nlohmann_json: \$(pkg-config --modversion nlohmann_json)
    kseq++: \$(pkg-config --modversion kseq++)
    bamtools: \$(pkg-config --modversion bamtools-1)
    eigen: 3.2
    boost: 1.84.0
END_VERSIONS
"""
}
