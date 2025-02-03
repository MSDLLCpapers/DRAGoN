
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process DownsampleBAM {
    label 'process_single'
    input:
        tuple val(meta), path(mergedbam), path(barcodes), path(demux_log)
    output:
        tuple val(meta), path("*_unmap.bam"), path(barcodes), optional: true, emit: bam
        tuple val(meta), path("${prefix}.downsample.txt"), emit: report
        tuple val(meta), path("*_discard.bam"), optional: true, emit: discard
        path 'versions.yml', emit: versions
    script:
        prefix = task.ext.prefix ?: meta.id
        args = task.ext.args ?: ''
"""
downsample_bam.sh \\
    -i ${mergedbam} \\
    -b ${barcodes} \\
    -s ${meta.id} \\
    -d ${demux_log} \\
    -c ${task.cpus} \\
    -p ${prefix} \\
    ${args}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    python: \$(python --version | sed -e "s/Python //g")
    numpy: \$(python -c 'import numpy; print(numpy.__version__)')
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: meta.id
"""
touch ${prefix}_unmap.bam
touch ${prefix}_discard.bam
if [ "${meta.bcidx}" != "noBC" ]; then
    i=0
    while IFS=\$'\t' read -r -a line; do
        if [ "\$i" -eq "${meta.bcidx}" ]; then
            break
        fi
        ((++i))
    done <"${barcodes}"
else
    i=-1
    line=(NNNNNNNNNNNNNNNN noBC -1 -1 True)
fi
echo -e "sample\treplicate\tindex\twell_name\tqc_pass_reads\tlimit\tseed\tfraction\tout\n${meta.id}\t${meta.id2}\t${meta.bcidx}\t\${line[1]}\t0\t0\t0\t0\t0" > ${prefix}.downsample.txt
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    python: \$(python --version | sed -e "s/Python //g")
    numpy: \$(python -c 'import numpy; print(numpy.__version__)')
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
END_VERSIONS
"""
}
