
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process Deduplicate {
    label 'process_single'
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
        if (true /*params.implementation == "C++" && task.attempt == 1*/) {
            program = 'deduplicate'
            // Eigen3 version is hard-coded because there's no CLI way to extract it
            versions = """
            "${task.process}":
                g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
                cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
                pkg-config: \$(pkg-config --version)
                zlib: \$(pkg-config --modversion zlib)
                bamtools: \$(pkg-config --modversion bamtools-1)
                eigen: 3.2
            """.stripIndent()
        } else {
            program = "collapse_umis.py"
            versions = """
            "${task.process}":
                python: \$(python --version | sed -e "s/Python //g")
                cutadapt: \$(cutadapt --version)
                biopython: \$(python -c 'import Bio; print(Bio.__version__)')
                pysam: \$(python -c 'import pysam; print(pysam.__version__)')
                numpy: \$(python -c 'import numpy; print(numpy.__version__)')
                pandas: \$(python -c 'import pandas; print(pandas.__version__)')
                scipy: \$(python -c 'import scipy; print(scipy.__version__)')
                xopen: \$(python -c 'import xopen; print(xopen.__version__)')
            """.stripIndent()
        }
"""
${program} "$annotatedBam" \\
    $barcodes \\
    $annofile \\
    "DRAGoN.out.${prefix}" \\
    --cpus ${task.cpus} \\
    --mmap-strategies $ambiguous \\
    ${args}
cat <<-END_VERSIONS > versions.yml
$versions
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
        if (true /*params.implementation == "C++" && task.attempt == 1*/) {
            program = 'deduplicate'
            // Eigen3 version is hard-coded because there's no CLI way to extract it
            versions = """
            "${task.process}":
                g++: \$(g++ --version | head -n1 | grep -oE '\\S+\$')
                cmake: \$(cmake --version | head -n1 | grep -oE '\\S+\$')
                pkg-config: \$(pkg-config --version)
                zlib: \$(pkg-config --modversion zlib)
                bamtools: \$(pkg-config --modversion bamtools-1)
                eigen: 3.2
            """.stripIndent()
        } else {
            program = "collapse_umis.py"
            versions = """
            "${task.process}":
                python: \$(python --version | sed -e "s/Python //g")
                cutadapt: \$(cutadapt --version)
                biopython: \$(python -c 'import Bio; print(Bio.__version__)')
                pysam: \$(python -c 'import pysam; print(pysam.__version__)')
                numpy: \$(python -c 'import numpy; print(numpy.__version__)')
                pandas: \$(python -c 'import pandas; print(pandas.__version__)')
                scipy: \$(python -c 'import scipy; print(scipy.__version__)')
                xopen: \$(python -c 'import xopen; print(xopen.__version__)')
            """.stripIndent()
        }
"""
touch ${annotatedBam.baseName}_MarkDups.bam
mkdir DRAGoN.out.${prefix}
touch DRAGoN.out.${prefix}/{matrix.mtx,features.tsv,barcodes.tsv}
cat <<-END_VERSIONS > versions.yml
$versions
END_VERSIONS
"""
}
