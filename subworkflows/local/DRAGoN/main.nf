
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {CountUnmappedReads} from "../../../modules/local/custom/count_unmapped_reads/main.nf"
include {NormalizeCountsWF} from "../../../subworkflows/local/NormalizeCounts/main.nf"
include {STAR_ALIGN} from '../../../modules/local/star/align/main'
include {FeatureCounts} from '../../../modules/local/subread/featurecounts/main'
include {SAMTOOLS_SORT as SortByName;
         SAMTOOLS_SORT as SortByCoord} from '../../../modules/local/samtools/sort/main'
include {Deduplicate} from '../../../modules/local/custom/deduplicate/main'
include {MergeDragonOutput} from '../../../modules/local/custom/merge_dragon/main'

workflow DRAGoN {
    take:
        decontam
        demux_log
        staridx
        annofile
        downsample_report
        val_ambiguous
        val_adjust_effective_lengths
        val_keep_bam

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()

        STAR_ALIGN (decontam, staridx, annofile, false)
        versions = versions.mix(STAR_ALIGN.out.versions.first())
        multiqc_files = multiqc_files.mix(STAR_ALIGN.out.log.map{ it -> it[1] })

        FeatureCounts (
            STAR_ALIGN.out.bam,
            annofile
        )
        versions = versions.mix(FeatureCounts.out.versions.first())
        multiqc_files = multiqc_files.mix(FeatureCounts.out.stats.map{ it -> it[1] })

        SortByName(FeatureCounts.out.bam)
        versions = versions.mix(SortByName.out.versions.first())

        Deduplicate (
            SortByName.out.bam,
            annofile,
            val_ambiguous
        )
        versions = versions.mix(Deduplicate.out.versions.first())

        if (val_adjust_effective_lengths || val_keep_bam) {
            SortByCoord(Deduplicate.out.bam)
            versions = versions.mix(SortByCoord.out.versions.first())
            bam_out = SortByCoord.out.bam
        } else {
            bam_out = Deduplicate.out.bam
        }
        CountUnmappedReads (
            STAR_ALIGN.out.unmapped
                .map { meta, bam -> [meta.id, meta, bam] }
                .groupTuple()
                .map { _id, meta, bam ->
                    def x = [meta, bam]
                        .transpose()
                        .sort{ it -> it[0].bcidx }
                        .transpose()
                    x[0] = x[0][0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)}
                    x
                }
        )
        versions = versions.mix(CountUnmappedReads.out.versions.first())

        merged_output = Deduplicate.out.counts
            .join(FeatureCounts.out.stats)
            .map{ meta, matrices, fcstats ->
                [ meta.id, meta, matrices, fcstats ]
            }
            .groupTuple()
            .join(demux_log.map{ meta, json ->
                [ meta.id, json ]
            })
            .join(CountUnmappedReads.out.report.map{ meta, report ->
                [ meta.id, report ]
            })
            .join(downsample_report.map{ meta, report ->
                [ meta.id, report ]
            })
            .map { _id, meta, matrices, fcstats, json, unmap, downsample ->
                def x = [meta, matrices, fcstats]
                    .transpose()
                    .sort{ it -> it[0].bcidx }
                    .transpose()
                x[0] = x[0][0].findAll{k, _v -> !['bcidx', 'bcname'].contains(k)}
                return x + [json, unmap, downsample]
            }
        MergeDragonOutput(
            merged_output,
            val_ambiguous
        )
        versions = versions.mix(MergeDragonOutput.out.versions.first())

        counts_matrix = MergeDragonOutput.out.merged
            .map{ meta, _report, _fcstats, _dedup, dense, _sparse ->
                [ meta, dense ]
            }
            .dump( tag: 'counts_matrix' )
        NormalizeCountsWF(
            Deduplicate.out.bam,
            counts_matrix,
            annofile,
            val_adjust_effective_lengths
        )
        versions = versions.mix(NormalizeCountsWF.out.versions)

    emit:
        marked_bam = bam_out.map{ meta, bam, _barcodes ->
            [meta, bam]
        }
        unmapped = STAR_ALIGN.out.unmapped
        counts_matrix
        normed_counts = NormalizeCountsWF.out.normed_counts
        versions
        multiqc_files
}
