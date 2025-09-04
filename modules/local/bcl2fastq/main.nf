// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include { usingProfile } from '../../../subworkflows/local/util_drugseq/main'

process BCL2FASTQ {
    label 'process_medium'
    container 'nf-core/bcl2fastq:2.20.0.422'

    input:
    path bcldir
    path sampleSheet

    output:
    path ("output/**_S[1-9]*_R?_00?.fastq.gz"), emit: fastq
    path ("output/**_S[1-9]*_I?_00?.fastq.gz"), optional: true, emit: fastq_idx
    path ("output/**Undetermined_S0*_R?_00?.fastq.gz"), optional: true, emit: undetermined
    path ("output/**Undetermined_S0*_I?_00?.fastq.gz"), optional: true, emit: undetermined_idx
    path ("output/Reports"), emit: reports
    path ("output/Stats"), emit: stats
    path ("InterOp/*.bin"), emit: interop
    path ("versions.yml"), emit: versions

    script:
    args = task.ext.args ?: ''
    nonIOcpus = [task.cpus - 2, 1].max()
    if (usingProfile('conda')) {
        error("process ${task.process} is not compatible with the 'conda' profile")
    }
    """
if [[ "${bcldir}" =~ .tar.gz\$ ]]; then
    tar xf "${bcldir}"
    bcldir="${bcldir.name.replaceAll(/\.tar\.gz$/, '')}"
else
    bcldir="${bcldir}"
fi
bcl2fastq -R \$bcldir -o output -r 1 -p ${nonIOcpus} -w 1 --sample-sheet ${sampleSheet} ${args}
cp -r \${bcldir}/InterOp .
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bcl2fastq: \$(bcl2fastq --version 2>&1 | grep '^bcl2fastq' | sed 's/^.+ v//')
END_VERSIONS
"""

    stub:
    args = task.ext.args ?: ''
    nonIOcpus = [task.cpus - 2, 1].max()
    """
# Get flowcell id
flowcell=\$(grep -Em1 '<Flowcell>\\w+</Flowcell>' ${bcldir}/RunInfo.xml | sed -r 's/<\\/?Flowcell>//')
mkdir -p output/\$flowcell Reports Stats InterOp

# Parse samplesheet to get output fastq paths
awk -F, -v basedir=output/\$flowcell -f - ${sampleSheet} <<END_AWK > paths.txt
BEGIN {k=1}
\\\$0 ~ /^\\\\[/ {
    isData = \\\$0 ~ /^\\[Data\\]/
    next
}
\\\$0 ~ /^\\\$/ {isData = false; next}
!isData {next}
!hasHeader {for (i = 1; i <= NF; i++) {header[\\\$i] = i} hasHeader = 1; next}
hasHeader {
    if ("Sample_Project" in header) {
        print basedir "/" \\\$header["Sample_Project"] "/" \\\$header["Sample_ID"] "_S" k "_L00" \$header["Lane"] "_R1_001.fastq.gz"
        print basedir "/" \\\$header["Sample_Project"] "/" \\\$header["Sample_ID"] "_S" k "_L00" \$header["Lane"] "_R2_001.fastq.gz"
    } else {
        print basedir "/" \\\$header["Sample_ID"] "_S" k "_L00" \\\$header["Lane"] "_R1_001.fastq.gz"
        print basedir "/" \\\$header["Sample_ID"] "_S" k "_L00" \\\$header["Lane"] "_R2_001.fastq.gz"
    }
    k++
}
END_AWK
while read path; do
    mkdir -p \$(dirname \$path)
    echo -ne '' | gzip -c > \$path
done <paths.txt

touch Reports/DUMMY Stats/DUMMY InterOp/DUMMY.bin
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bcl2fastq: \$(bcl2fastq --version 2>&1 | grep '^bcl2fastq' | sed 's/^.+ v//')
END_VERSIONS
"""
}
