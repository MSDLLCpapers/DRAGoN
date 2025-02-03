
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {SAMTOOLS_STATS} from '../../../modules/local/samtools/stats/main'
include {GenNormalizedCountsTable} from '../../../modules/local/custom/normalize_counts/main'

workflow NormalizeCountsWF {
    take:
        bams
        raw_counts
        annofile
        val_adjust_effective_lengths

    main:
        versions = Channel.empty()
        if (val_adjust_effective_lengths) {
            SAMTOOLS_STATS(bams)
            fragment_lengths = SAMTOOLS_STATS.out.stats
                .map{ meta, stats ->
                    def length = stats
                        .readLines()
                        .find{line -> line.startsWith('SN\taverage length:')}
                        .strip()
                        .tokenize("\t")
                        .getAt(2)
                    return [ meta.id, meta, length as Integer ]
                }
                .groupTuple()
            versions = versions.mix(SAMTOOLS_STATS.out.versions.first())
        } else {
            fragment_lengths = raw_counts.map{ meta, _counts ->
                [meta.id, [], []]
            }
        }
        ch_norm_in = raw_counts
            .map {meta, counts ->
                [ meta.id, meta, counts ]
            }
            .combine(fragment_lengths, by: 0)
            .flatMap { _id, meta, counts, meta2, length ->
                def flens = [meta2, length]
                    .transpose()
                    .sort{ it -> it[0].bcidx }
                    .transpose()
                return counts.collect{ ct ->
                    [ meta + [ bcidxs: flens[0].collect{ it -> it.bcname }], ct, flens[1] ]
                }
            }
            .dump(tag: 'ch_norm_in')
        GenNormalizedCountsTable(
            ch_norm_in.combine(['tpm', 'fpkm']),
            annofile.collect(),
        )
        versions = versions.mix(GenNormalizedCountsTable.out.versions.first())

    emit:
        normed_counts = GenNormalizedCountsTable.out.counts
        versions
}
