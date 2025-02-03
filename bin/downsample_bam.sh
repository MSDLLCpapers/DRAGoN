#!/usr/bin/env bash

# Written by Scott Norton
# 2024-06-27

# Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
# This file is part of DRAGoN.
#
# This source code is licensed under the MIT License found in the
# LICENSE file in the root directory of this source tree.


set -euo pipefail

which jq &>/dev/null || { echo -e "jq: no such file or directory"; exit 127; }
which samtools &>/dev/null || { echo -e "samtools: no such file or directory"; exit 127; }
which python &>/dev/null || { echo -e "python: no such file or directory"; exit 127; }
python -c "import pandas as pd" &>/dev/null || { echo -e "python: unable to import pandas: no such package"; exit 127; }

INT64_MAX=9223372036854775807

inbam=
barcodes=
samplename=
prefix=
inlog=
threshold=false
downsample=false
cpus=1
bcidx=noBC
rngseed=
usage="usage: $(basename "$0") [OPTS]

  options:
    -i <BAM>        Input BAM
    -b <TSV>        Barcodes file with 5 columns: seq, name, x, y, mask
    -n <STR>        Barcode index (def. $bcidx)
    -s <STR>        Sample name
    -p <STR>        Output prefix
    -d <JSON>       Demux JSON
    -t [INT|FLT]    Downsampling target (int: # of reads; float: fraction of well size; no arg: dynamic; default: no downsampling)
    -T [INT|FLT]    Downsampling threshold (int: # of reads; float: fraction of experiment size; no arg: dynamic; default: no downsampling)
    -c <INT>        Number of extra CPUs (def. $cpus)
    -k <INT>        RNG seed"

while getopts ":i:b:n:s:p:tTd:c:k:h" opt; do
    case $opt in
        h )
            echo -e "$usage"
            exit 0 ;;
        i )
            inbam="$OPTARG" ;;
        b )
            barcodes="$OPTARG" ;;
        n )
            bcidx="$OPTARG" ;;
        s )
            samplename="$OPTARG" ;;
        p )
            prefix="$OPTARG" ;;
        t )
            nextopt=${*:$OPTIND:1}
            if [[ -n $nextopt && $nextopt != -* ]]; then
                ((++OPTIND))
                downsample=$nextopt
            else
                downsample=true
            fi ;;
        T )
            nextopt=${*:$OPTIND:1}
            if [[ -n $nextopt && $nextopt != -* ]]; then
                ((++OPTIND))
                threshold=$nextopt
            else
                threshold=true
            fi ;;
        c )
            cpus="$OPTARG" ;;
        d )
            inlog="$OPTARG" ;;
        k )
            rngseed="$OPTARG" ;;
        : )
            echo -e "Invalid option: -$OPTARG requires an argument\n$usage" 1>&2
            exit 255 ;;
        \? )
            echo -e "Error: Invalid option\n$usage" 1>&2
            exit 255 ;;
    esac
done
shift $((OPTIND-1))
set +u
if [ -n "$1" ]; then
    echo -e "Invalid args: $*\n$usage" 1>&2
    exit 255
fi
set -u

# check for missing or conflicting args
missing() {
    echo -e "Missing or invalid required argument: $1\n$usage" 1>&2
    exit 255
}
invalid() {
    echo -e "Invalid optional argument: $1\n$usage" 1>&2
    exit 255
}
args_error() {
    echo -e "Erroneous arguments: $1\n$usage" 1>&2
    exit 255
}
if [ -z "$inbam" ] || ! [ -f "$inbam" ]; then
    missing "-i <BAM>"
fi
if [ -z "$barcodes" ] || ! [ -f "$barcodes" ]; then
    missing "-b <TSV>"
fi
if [ -z "$samplename" ]; then
    misisng "-s <STR>"
fi
if [ -z "$prefix" ]; then
    misisng "-p <STR>"
fi
if [ -z "$inlog" ]; then
    missing "-d <JSON>"
fi
int_re='^[0-9]+$'
if [[ ! "$cpus" =~ $int_re ]]; then
    invalid "-c <INT>"
fi
downsample_re='^(true|false|(0?\.)?[0-9]+)$'
if [[ ! "$downsample" =~ $downsample_re ]]; then
    invalid "-t [INT|FLT]"
fi
if [[ ! "$threshold" =~ $downsample_re ]]; then
    invalid "-T [INT|FLT]"
fi
if [ "$downsample" = "false" ] && [ "$threshold" != "false" ]; then
    downsample="$threshold"
elif [ "$downsample" != "false" ] && [ "$threshold" = "false" ]; then
    threshold="$downsample"
fi
if [ -n "$rngseed" ]; then
    if [[ ! "$rngseed" =~ $int_re ]]; then
        invalid "-k <INT>"
    fi
    RANDOM=$rngseed
fi

experiment_size=$(jq ".read_count" "$inlog")

calc_dynthresh() {
    echo "Computing dynamic threshold" 1>&2
    iqrvec=()
    num_qc_pass=0
    while IFS=$'\t' read -r -a line; do
        if [ "${#line[*]}" -lt 5 ] || [ "${line[4]}" != "True" ]; then
            addend=$(jq ".fail_counts.PASS.\"${line[1]}\"" "$inlog")
            ((num_qc_pass+=addend))
            iqrvec+=("$addend")
        fi
    done < "$barcodes"
    if [ "$num_qc_pass" -eq 0 ]; then
        echo "ERROR: No reads from unmasked wells passing QC" 1>&2
        exit 1
    fi
    python - <<EOF
import pandas as pd
data = pd.Series("${iqrvec[*]}".split(' ')).astype(int)
q25, q75 = data.quantile((0.25, 0.75))
print(min(2.75 * q75 - 1.75 * q25, ${experiment_size} * 0.15))
EOF
}

dynthresh=
case $threshold in
    false )
        limit=$INT64_MAX ;;
    true )
        dynthresh=$(calc_dynthresh)
        limit=$(python -c "print(round($dynthresh))") ;;
    ?.* )
        limit=$(python -c "print(round(${experiment_size} * ${threshold}))") ;;
    * )
        limit=$threshold ;;
esac
if [ "$limit" -eq "$INT64_MAX" ]; then
    lim_out=inf
else
    lim_out="$limit"
fi
echo -e "sample\treplicate\tindex\twell_name\tqc_pass_reads\tlimit\tseed\tfraction\tout" >"${prefix}.downsample.txt"
verbose=1
case $downsample in
    false )
        verbose=0
        ds_msg="so it is disabled" ;;
    true )
        if [ -z "$dynthresh" ]; then
            dynthresh=$(calc_dynthresh)
        fi
        ds_msg="to a target level of $dynthresh reads across lanes" ;;
    ?.* )
        ds_msg="to a target fixed ratio of $downsample" ;;
    * )
        ds_msg="to a target level of $downsample reads across lanes" ;;
esac
echo "Downsampling with a threshold of $lim_out $ds_msg"
if [ "$bcidx" != "noBC" ]; then
    i=0
    while IFS=$'\t' read -r -a line; do
        if [ "$i" -eq "$bcidx" ]; then
            break
        fi
        ((++i))
    done <"$barcodes"
else
    i=-1
    line=(NNNNNNNNNNNNNNNN noBC -1 -1 True)
fi
keepbam="${prefix}_unmap.bam"
dropbam="${prefix}_discard.bam"
seed=$RANDOM
num_in=$(jq ".fail_counts.PASS.\"${line[1]}\"" "$inlog")
cur_lane_num_in=$(samtools view -@ "$cpus" -F QCFAIL -c "$inbam")
special=
if [ "${#line[*]}" -ge 5 ] && [ "${line[4]}" = "True" ]; then
    fraction=0
    special=" because it is softmasked"
elif [ "$num_in" -le "$limit" ]; then
    fraction=1
    special=" because it is below the downsampling threshold"
else
    case $downsample in
        false )
            fraction=1 ;;
        true )
            fraction=$(python -c "print(min(1, $limit / $num_in))") ;;
        ?.* )
            fraction=$downsample ;;
        * )
            fraction=$(python -c "print(min(1, $downsample / $num_in))") ;;
    esac
fi

if [[ "$fraction" =~ ^1(\.[0-9]+)?$ ]]; then
    if [ $verbose -eq 1 ]; then
        echo "Well ${line[1]} is not downsampled$special"
    fi
    samtools view -@ "$cpus" -F QCFAIL -U "$dropbam" -o "$keepbam" "$inbam"
    num=$cur_lane_num_in
elif [[ "$fraction" =~ ^0(\.0)?$ ]]; then
    if [ $verbose -eq 1 ]; then
        echo "Well ${line[1]} is downsampled to 0$special"
    fi
    cp "$inbam" "$dropbam"
    num=0
else
    samtools view -@ "$cpus" -F QCFAIL --subsample "$fraction" --subsample-seed "$seed" -U "$dropbam" -o "$keepbam" "$inbam"
    num=$(samtools view -c -@ "$cpus" "$keepbam")
    if [ $verbose -eq 1 ]; then
        echo "Well ${line[1]} was downsampled to a ratio of ${fraction} ($(cat tmp.txt) --> ${num})"
    fi
fi
echo -e "${samplename}\t${prefix}\t${i}\t${line[1]}\t$cur_lane_num_in\t${lim_out}\t${seed}\t${fraction}\t${num}" >>"${prefix}.downsample.txt"
