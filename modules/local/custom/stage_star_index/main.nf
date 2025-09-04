process StageStarIndex {
    executor 'local'
    label 'process_single'
    stageInMode 'copy'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:24.04'
        : 'biocontainers/ubuntu:24.04'}"

    input:
    path files

    output:
    path files, emit: files
    path 'versions.yml', emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: 24.04
    END_VERSIONS
    """
}
