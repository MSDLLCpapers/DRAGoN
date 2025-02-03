
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {SAMTOOLS_MERGE} from '../../../modules/local/samtools/merge/main'

workflow Postprocessing {
    take:
        annobam
        contamBam
        unmappedBam
        discardedBam
        val_keep_bam

    main:
        versions = Channel.empty()
        ch_bam = Channel.empty()
        ch_bai = Channel.empty()
        if (val_keep_bam) {
            Channel.empty()
                | mix(annobam, contamBam, unmappedBam, discardedBam)
                | transpose
                | map{ meta, bam ->
                    [ meta.id, meta, bam ]
                }
                | groupTuple
                | map { _id, meta, bam ->
                    [
                        meta[0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)},
                        bam.sort{b -> b.name}
                    ]
                }
                | SAMTOOLS_MERGE
            versions = versions.mix(SAMTOOLS_MERGE.out.versions.first())
            ch_bam = SAMTOOLS_MERGE.out.bam
            ch_bai = SAMTOOLS_MERGE.out.bai
        }

    emit:
        ch_bam
        ch_bai
        versions
}
