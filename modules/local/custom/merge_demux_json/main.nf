
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process MergeDemuxJson {
    label 'process_single'
    input:
        tuple val(meta), path(json, name: "input??/*"), path(barcodes)
        path plate_layout
        val data_quality_errors
    output:
        tuple val(meta), path("demux.json"), stdout, emit: merged
        path "*.png", emit: pngs
        path 'versions.yml', emit: versions
    script:
        args = task.ext.args ?: ""
"""
merge_demux_json.py \\
    --handle-outlier-wells ${data_quality_errors} \\
    --layout "${plate_layout}" \\
    --barcodes "${barcodes}" \\
    --nobc-resolve-threads "${task.cpus}" \\
    ${args} \\
    ${json}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    seaborn: \$(python -c 'import seaborn; print(seaborn.__version__)')
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    adjustText: \$(python -c 'import adjustText; print(adjustText.__version__)')
END_VERSIONS
"""
    stub:
"""
touch demux.json violin.png heatmap.png
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    seaborn: \$(python -c 'import seaborn; print(seaborn.__version__)')
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    adjustText: \$(python -c 'import adjustText; print(adjustText.__version__)')
END_VERSIONS
"""
}
