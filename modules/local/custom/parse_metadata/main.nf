// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process ParseMetadataXLSX {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-9f99278ff296588175d0d58987544943664a125d:ba896f1b212b989d446396b179f05708fe657002-0'
        : 'biocontainers/mulled-v2-9f99278ff296588175d0d58987544943664a125d:ba896f1b212b989d446396b179f05708fe657002-0'}"

    input:
    path metadata
    path samplesheet

    output:
    path "*_barcodes.txt", emit: barcodes
    path samplesheet, emit: sampleSheet
    path 'versions.yml', emit: versions

    script:
    sampleSheetArg = samplesheet ? "-s ${samplesheet}" : ''
    """
metadata_to_run_specs.py ${metadata} ${sampleSheetArg}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    openpyxl: \$(python -c 'import openpyxl; print(openpyxl.__version__)')
END_VERSIONS
"""

    stub:
    sampleSheetArg = samplesheet ? "-s ${samplesheet}" : ''
    """
metadata_to_run_specs.py ${metadata} ${sampleSheetArg}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    openpyxl: \$(python -c 'import openpyxl; print(openpyxl.__version__)')
END_VERSIONS
"""
}
