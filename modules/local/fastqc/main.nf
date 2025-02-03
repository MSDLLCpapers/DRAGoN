
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process FastQC {
    label 'process_medium'
    input:
        tuple val(meta), path(inputs, arity: '1..*')
        val mapped
    output:
        path "${prefix}", emit: zip
        path 'versions.yml', emit: versions
    script:
        prefix = task.ext.prefix ?: "${meta.stage?:'untrimmed'}.${meta.id2?:''}"
        filetype = meta.filetype ?: (inputs[0].extension == 'gz' ? inputs[0].baseName.tokenize(/\./)[-1] : inputs[0].extension)
        assert ["bam", "sam", "fastq"].contains(filetype)
        if (filetype != "fastq" && mapped) {
            filetype = "${filetype}_mapped"
        }
"""
mkdir -p ${prefix}
fastqc -o ${prefix} --extract -f ${filetype} -t ${task.cpus} --memory ${[task.memory.mega, 10000].min()} \
    --casava $inputs
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastqc: \$(fastqc --version | sed 's/FastQC v//g')
END_VERSIONS
"""
    stub:
        prefix = task.ext.prefix ?: "${meta.stage?:'untrimmed'}.${meta.id2?:''}"
"""
mkdir -p ${prefix}
touch ${inputs.collect{name -> prefix + '/' + name + '.zip'}.join(' ')}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastqc: \$(fastqc --version | sed 's/FastQC v//g')
END_VERSIONS
"""
}
