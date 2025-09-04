// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

// STAR_ALIGN_SE
process STAR_ALIGN {
    label 'process_high'
    label 'process_high_memory'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'}"

    input:
    tuple val(meta), path(unmapSam, name: "input?/*", arity: '1..*'), path(barcodes)
    path idx, name: 'idx/*'
    path annofile
    val isSolo

    output:
    tuple val(meta), path("*Aligned.out.bam"), path(barcodes), optional: true, emit: bam
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam"), optional: true, emit: bam_sorted
    tuple val(meta), path("*Solo.out"), path(barcodes), optional: true, emit: solo
    tuple val(meta), path("*Unmapped.out.bam"), optional: true, emit: unmapped
    tuple val(meta), path("*Log.final.out"), emit: log
    path 'versions.yml', emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    args = task.ext.args ?: ''
    gtfFileArg = annofile ? "--sjdbGTFfile ${annofile}" : ""
    """
if [ ${unmapSam.size()} -eq 1 ]; then
    ln -s ${unmapSam} ${prefix}.in.bam
else
    samtools cat -@ ${task.cpus} -o ${prefix}.in.bam ${unmapSam}
fi
STAR \\
    --runThreadN "${task.cpus}" \\
    --readFilesIn "${prefix}.in.bam" \\
    --outFileNamePrefix "${prefix}_" \\
    --genomeDir idx \\
    ${gtfFileArg} \\
    ${args}

if [ -f ${prefix}_Unmapped.out.mate1 ]; then
    (samtools view -H "${prefix}.in.bam"; samtools import -T '*' -0 "${prefix}_Unmapped.out.mate1") | samtools view -b -@ ${task.cpus} -o "${prefix}_Unmapped.out.bam"
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
END_VERSIONS
"""

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.bcidx}"
    """
touch ${prefix}_Aligned.out.bam
touch ${prefix}_Aligned.sortedByCoord.out.bam
touch ${prefix}_Log.final.out
touch ${prefix}_Unmapped.out.bam
mkdir ${prefix}_Solo.out
touch ${prefix}_Solo.out/matrix.mtx
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    star: \$(STAR --version | sed -e "s/STAR_//g")
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
END_VERSIONS
"""
}
