include { BCFTOOLS_STATS                  } from '../../../modules/nf-core/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL  } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_LDEPTH  } from '../../../modules/nf-core/vcftools/main'

workflow VCF_QC_BCFTOOLS_VCFTOOLS {
    take:
        vcf
        tbi
        target_bed

    main:

    ch_versions = Channel.empty()
    
    vcf_maybe_tbi = tbi ? : vcf.join(tbi) : vcf
    BCFTOOLS_STATS(vcf_maybe_tbi target_bed, [], [])
    VCFTOOLS_TSTV_COUNT(vcf, target_bed, [])
    VCFTOOLS_TSTV_QUAL(vcf, target_bed, [])
    VCFTOOLS_SUMMARY(vcf, target_bed, [])
    VCFTOOLS_LDEPTH(vcf, target_bed, [])

    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_TSTV_COUNT.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_counts    = VCFTOOLS_TSTV_COUNT.out.tstv_count
    vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
    vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary
    vcftools_ldepth         = VCFTOOLS_LDEPTH.out.ldepth
    versions                = ch_versions
}
