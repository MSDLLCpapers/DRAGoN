// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process PIGZ {
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pigz:2.8'
        : 'biocontainers/pigz:2.8'}"
    stageInMode 'copy'  // copy instead of symlink, for pigz requirements

    input:
    tuple val(meta), path(file_in)

    output:
    tuple val(meta), path("${prefix}${suffix}"), emit: result
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: file_in.baseName
    suffix = task.ext.suffix ?: ''
    args = task.ext.args ?: ''
    if (args =~ /-\w*d\w*/ && file_in.extension != 'gz') {
        error("expected a gzip file, got ${file_in}")
    }
    if (file_in.name == "${prefix}.${suffix}") {
        error("input and output filenames are the same")
    }
    pipe = args =~ /-\w*c\w*/ ? "> ${prefix}${suffix}" : ''
    """
        pigz \\
            ${args} \\
            -p ${task.cpus} \\
            ${file_in} \\
            ${pipe}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pigz: \$(pigz --version | cut -d' ' -f2)
        END_VERSIONS
        """

    stub:
    prefix = task.ext.prefix ?: file_in.baseName
    suffix = task.ext.suffix ?: ''
    if (file_in.name == "${prefix}.${suffix}") {
        error("input and output filenames are the same")
    }
    """
        ln -s file_in ${prefix}.${suffix}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pigz: \$(pigz --version | cut -d' ' -f2)
        END_VERSIONS
        """
}
