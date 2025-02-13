
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

process GenerateReport {
    label 'process_single'
    input:
        tuple val(meta), path(solo_outs), path(barcodes), path(demux_logs), path(unmapped), path(downsample)
        val ambiguous
    output:
        tuple val(meta), path("STARsoloReport.txt"), path("counts*.tsv", arity: '1..5'), path("Solo.out"), emit: report
        path 'versions.yml', emit: versions
    script:
"""
starsolo_report.py --demux ${demux_logs} --solo-outs ${solo_outs} \
    --barcodes ${barcodes} --ambiguous ${ambiguous} \
    --unmapped-report ${unmapped} --downsample-report $downsample
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    scipy: \$(python -c 'import scipy; print(scipy.__version__)')
END_VERSIONS
"""
    stub:
"""
touch STARsoloReport.txt counts.tsv
mkdir -p Solo.out
touch Solo.out/{Gene.stats,matrix.mtx,features.tsv,barcodes.tsv}
for amb in ${ambiguous.tokenize(/ /).findAll{amb -> amb != 'Unique'}.join(' ')}; do
    touch counts-\$amb.tsv Solo.out/UniqueAndMult-\$amb.mtx
done
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed -e "s/Python //g")
    pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    scipy: \$(python -c 'import scipy; print(scipy.__version__)')
END_VERSIONS
"""
}
