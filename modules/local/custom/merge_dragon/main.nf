// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process MergeDragonOutput {
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-f5bcc49db2e8d41c3472bc591a910399d68ee50e:0f862c589fc95bcef32fecaa39ca79a7f5b2abf1-0'
        : 'biocontainers/mulled-v2-f5bcc49db2e8d41c3472bc591a910399d68ee50e:0f862c589fc95bcef32fecaa39ca79a7f5b2abf1-0'}"

    input:
    tuple val(meta), path(counts), path(feature_stats), path(dmux_stats), path(unmapped), path(downsample)
    val ambiguous

    output:
    tuple val(meta), path("DRAGoNreport.txt"), path("featureCounts.summary"), path("dedup_stats.tsv"), path("counts*.tsv", arity: '1..5'), path("DRAGoN.out"), emit: merged
    path 'versions.yml', emit: versions

    script:
    """
dragon_report.py --matrices ${counts} --dmux ${dmux_stats} \
    --features ${feature_stats} --ambiguous ${ambiguous} \
    --unmapped-report ${unmapped} --downsample-report ${downsample}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    scipy: \$(python -c 'import scipy; print(scipy.__version__)')
END_VERSIONS
"""

    stub:
    """
touch DRAGoNreport.txt featureCounts.summary dedup_stats.tsv counts.tsv
mkdir DRAGoN.out
touch DRAGoN.out/{matrix.mtx,features.tsv.gz,barcodes.tsv.gz}
for amb in ${ambiguous.tokenize(/ /).findAll { amb -> amb != 'Unique' }.join(' ')}; do
    touch counts-\$amb.tsv DRAGoN.out/UniqueAndMult-\$amb.mtx
done
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    scipy: \$(python -c 'import scipy; print(scipy.__version__)')
END_VERSIONS
"""
}
