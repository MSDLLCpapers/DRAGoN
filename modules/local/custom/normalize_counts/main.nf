
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process GenNormalizedCountsTable {
    label 'process_single'
    input:
        tuple val(meta), path(raw_counts), val(fraglen), val(kind)
        path annofile
    output:
        path "${prefix}.tsv", emit: counts
        path 'versions.yml', emit: versions
    script:
        prefix = task.ext.prefix ?: "${meta.id}.${raw_counts.baseName}.${kind}"
        fraglen_arg = meta.bcidxs && fraglen ? "-n ${meta.bcidxs.join(" ")} -f ${fraglen.join(" ")}" : ""
"""
counts_to_tpm.py $kind -c $raw_counts -a $annofile $fraglen_arg -o ${prefix}.tsv
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: "${meta.id}.${raw_counts.baseName}.${kind}"
"""
touch ${prefix}.tsv
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
END_VERSIONS
"""
}
