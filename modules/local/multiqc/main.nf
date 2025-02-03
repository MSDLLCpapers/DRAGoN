
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process MultiQC {
    label 'process_single'
    input:
        path inputs, stageAs: "?/*"
        path "multiqc_config.yaml"
    output:
        path "multiqc", emit: result
        path 'versions.yml', emit: versions
    script:
"""
multiqc -o multiqc .
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    multiqc: \$(multiqc --version | grep -oE '\\S+\$')
END_VERSIONS
"""
    stub:
"""
mkdir -p multiqc
touch multiqc/DUMMY
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    multiqc: \$(multiqc --version | grep -oE '\\S+\$')
END_VERSIONS
"""
}
