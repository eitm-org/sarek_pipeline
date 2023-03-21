//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_FIXVCFHEADER                 as FIX_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'
include { GATK4_FIXVCFHEADER                 as FIX_NORMAL_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'
include { GATK4_FIXVCFHEADER                 as FIX_TUMOR_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'

include { GATK4_MERGEVCFS                     as MERGE_VCFS_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_NORMAL_VCFS_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_TUMOR_VCFS_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'

include { CLAIRS                          as CLAIRS_PAIRED               } from '../../../modules/local/clairs'

workflow BAM_VARIANT_CALLING_SOMATIC_CLAIRS {
    take:
    input                     // channel: [ val(meta), [ input ], [ input_index ], [which_norm] ]
    fasta                     // channel: /path/to/reference/fasta
    fai                       // channel: /path/to/reference/fasta/index
    dict                      // channel: /path/to/reference/fasta/dictionary
    germline_resource         // channel: /path/to/germline/resource
    germline_resource_tbi     // channel: /path/to/germline/index
    vcf_header                // channel: /path/to/vcf_header
    normal_vcf                // channel: /path/to/normal_germline_vcf

    main:
    ch_versions = Channel.empty()

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //
    CLAIRS_PAIRED(
        input,
        fasta,
        fai,
        dict,
        normal_vcf
    )

    // Merge somatic VCF
    // Figure out if using intervals or no_intervals
    CLAIRS_PAIRED.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_vcf_branch }

    CLAIRS_PAIRED.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_tbi_branch }
    FIX_VCFHEADER_CLAIRS(
        clairs_vcf_branch.intervals
        .map{ meta, vcf -> [meta, vcf]},
        vcf_header
    )
    FIX_VCFHEADER_CLAIRS.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_fixed_vcf_branch }

    //Only when using intervals
    MERGE_VCFS_CLAIRS(
        clairs_fixed_vcf_branch.intervals
        .map{ meta, vcf ->

            new_meta = [
                        id:meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf]
        }.groupTuple(),
        dict
    )

    clairs_vcf = Channel.empty().mix(
        MERGE_VCFS_CLAIRS.out.vcf,
        clairs_vcf_branch.no_intervals)

    clairs_tbi = Channel.empty().mix(
        MERGE_VCFS_CLAIRS.out.tbi,
        clairs_tbi_branch.no_intervals)

    // Merge tumor germline VCF
    // Figure out if using intervals or no_intervals
    CLAIRS_PAIRED.out.vcf_germline_tumor.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_vcf_germline_tumor_branch }

    //Only when using intervals
    MERGE_TUMOR_VCFS_CLAIRS(
        clairs_vcf_germline_tumor_branch.intervals
        .map{ meta, vcf_tumor ->

            new_meta = [
                        id:meta.tumor_id + "_germline",
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf_tumor]
        }.groupTuple(),
        dict
    )

    clairs_vcf_germline_tumor = Channel.empty().mix(
        MERGE_TUMOR_VCFS_CLAIRS.out.vcf,
        clairs_vcf_germline_tumor_branch.no_intervals)



    // Merge normal germline VCF
    // Figure out if using intervals or no_intervals
    CLAIRS_PAIRED.out.vcf_germline_normal.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_vcf_germline_normal_branch }

    //Only when using intervals
    MERGE_NORMAL_VCFS_CLAIRS(
        clairs_vcf_germline_normal_branch.intervals
        .map{ meta, vcf_normal ->

            new_meta = [
                        id:meta.normal_id + "_germline",
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]

            [groupKey(new_meta, meta.num_intervals), vcf_normal]
        }.groupTuple(),
        dict
    )
    clairs_vcf_germline_normal = Channel.empty().mix(
        MERGE_NORMAL_VCFS_CLAIRS.out.vcf,
        clairs_vcf_germline_normal_branch.no_intervals)

    

    ch_versions = ch_versions.mix(MERGE_VCFS_CLAIRS.out.versions)
    ch_versions = ch_versions.mix(CLAIRS_PAIRED.out.versions)

    emit:
    clairs_vcf            = clairs_vcf                                     // channel: [ val(meta), [ vcf ] ]
    clairs_vcf_germline_normal = clairs_vcf_germline_normal                // channel: [ val(meta), [ vcf ] ]
    clairs_vcf_germline_tumor = clairs_vcf_germline_tumor                  // channel: [ val(meta), [ vcf ] ]
    versions               = ch_versions                                   // channel: [ versions.yml ]
}
