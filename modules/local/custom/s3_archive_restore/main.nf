
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process MaybeRestoreFromS3 {
    label 'process_high'
    label 'boto3'

    input:
        tuple val(meta), val(fastqfiles)
    output:
        tuple val(meta), val(fastqfiles), emit: fastq
        path 'versions.yml', emit: versions
    script:
"""
if [ ${!fastqfiles.flatten().isEmpty()} = true ]; then
    s3_restore_fastq.py -t ${task.cpus} ${fastqfiles.flatten().join(" ")}
fi
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    aioboto3: \$(python -c 'import aioboto3; print(aioboto3.__version__)')
END_VERSIONS
"""
    stub:
"""
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    aioboto3: \$(python -c 'import aioboto3; print(aioboto3.__version__)')
END_VERSIONS
"""
}
