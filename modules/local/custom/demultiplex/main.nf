
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

// Demultiplex the reads. STARsolo has its own
// demultiplexer but it's limited to 0 or 1 mismatches.
// cutadapt also has its own demultiplexer but it does not seem to work properly
// This approach accepts 2 or more mismatches and also
// splits the input into separate files for input.
process Demultiplex {
    label 'process_single'
    input:
        tuple val(meta), path(fastq), path(barcodes)
    output:
        tuple val(meta), path("${prefix}.bam"), path(barcodes), emit: bam
        tuple val(meta), path("*.demux.json"), path(barcodes), emit: log
        tuple val(meta), path("${prefix}_qcFail.bam"), optional:true, emit: qcfail
        path 'versions.yml', emit: versions
    script:
        prefix = "${task.ext.prefix ?: meta.id}"
        args = task.ext.args ?: ''
        if (true /*params.implementation == "C++" && task.attempt == 1*/) {
            program = 'demultiplex'
            // BOOST and Eigen3 versions are hard-coded because there's no CLI way to extract them
            versions = """
            "${task.process}":
                g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
                cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
                pkg-config: \$(pkg-config --version)
                zlib: \$(pkg-config --modversion zlib)
                kseq++: \$(pkg-config --modversion kseq++)
                bamtools: \$(pkg-config --modversion bamtools-1)
                eigen: 3.2
                boost: 1.84.0
            """.stripIndent()
        } else {
            program = "demultiplex.py"
            versions = """
            "${task.process}":
                python: \$(python --version | sed -e "s/Python //g")
                cutadapt: \$(cutadapt --version)
                pysam: \$(python -c 'import pysam; print(pysam.__version__)')
                regex: \$(python -c 'import regex; print(regex.__version__)')
                pandas: \$(python -c 'import pandas; print(pandas.__version__)')
            """.stripIndent()
        }
"""
set -ueo pipefail
${program} \\
    $fastq \\
    "$barcodes" \\
    "${prefix}" \\
    --cpus ${task.cpus} \\
    ${args}
cat <<-END_VERSIONS > versions.yml
$versions
END_VERSIONS
"""
    stub:
        prefix = "${task.ext.prefix ?: meta.id}"
        if (true /*params.implementation == "C++" && task.attempt == 1*/) {
            program = 'demultiplex'
            versions = """
            "${task.process}":
                g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
                cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
                pkg-config: \$(pkg-config --version)
                zlib: \$(pkg-config --modversion zlib)
                kseq++: \$(pkg-config --modversion kseq++)
                bamtools: \$(pkg-config --modversion bamtools-1)
                eigen: 3.2
                boost: 1.84.0
            """.stripIndent()
        } else {
            program = "demultiplex.py"
            versions = """
            "${task.process}":
                python: \$(python --version | sed -e "s/Python //g")
                cutadapt: \$(cutadapt --version)
                pysam: \$(python -c 'import pysam; print(pysam.__version__)')
                regex: \$(python -c 'import regex; print(regex.__version__)')
                pandas: \$(python -c 'import pandas; print(pandas.__version__)')
            """.stripIndent()
        }
"""
touch ${prefix}.bam ${prefix}_qcFail.bam ${prefix}.demux.json
cat <<-END_VERSIONS > versions.yml
$versions
END_VERSIONS
"""
}
