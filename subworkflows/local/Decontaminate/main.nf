
// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

include {STAR_ALIGN} from '../../../modules/local/star/align/main'

workflow Decontaminate {
    take:
        reads_in
        val_contaminant

    main:
        versions = Channel.empty()
        multiqc_files = Channel.empty()

        if (val_contaminant != null) {
            log.warn "Contaminant module is deprecated"
            val_contaminant = Channel.value(file(val_contaminant, type: 'dir', checkIfExists: true).listFiles() as List)
            STAR_ALIGN(
                reads_in,
                val_contaminant,
                [],
                false
            )
            versions = versions.mix(STAR_ALIGN.out.versions.first())
            multiqc_files = multiqc_files.mix(STAR_ALIGN.out.log)
        }

    emit:
        decontaminated = val_contaminant ? STAR_ALIGN.out.unmapped : reads_in
        contaminated = val_contaminant ? STAR_ALIGN.out.bam : Channel.empty()
        multiqc_files
        versions
}
