// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process GenNormalizedCountsTable {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-f5bcc49db2e8d41c3472bc591a910399d68ee50e:0f862c589fc95bcef32fecaa39ca79a7f5b2abf1-0'
        : 'biocontainers/mulled-v2-f5bcc49db2e8d41c3472bc591a910399d68ee50e:0f862c589fc95bcef32fecaa39ca79a7f5b2abf1-0'}"

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
counts_to_tpm.py ${kind} -c ${raw_counts} -a ${annofile} ${fraglen_arg} -o ${prefix}.tsv
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
