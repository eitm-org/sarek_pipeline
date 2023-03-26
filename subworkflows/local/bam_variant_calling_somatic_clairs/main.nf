//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_FIXVCFHEADER                 as FIX_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'
include { GATK4_FIXVCFHEADER                 as FIX_NORMAL_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'
include { GATK4_FIXVCFHEADER                 as FIX_TUMOR_VCFHEADER_CLAIRS               } from '../../../modules/local/gatk_fixvcfheader'

include { GATK4_MERGEVCFS                     as MERGE_VCFS_PAIRED_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_VCFS_TUMOR_GERMLINE_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_VCFS_TUMOR_PILEUP_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_VCFS_NORMAL_GERMLINE_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MERGEVCFS                     as MERGE_VCFS_NORMAL_PILEUP_CLAIRS               } from '../../../modules/nf-core/gatk4/mergevcfs'


include { CLAIRS                          as CLAIRS_PAIRED_FIRST               } from '../../../modules/local/clairs'
include { CLAIRS                          as CLAIRS_PAIRED_REST               } from '../../../modules/local/clairs'

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
    input.branch{
        first: it[0].normal_id[-1] == '0'
        rest:  it[0].normal_id[-1] != '0'
    }.set{input_branch}
    CLAIRS_PAIRED_FIRST(
        input_branch.first.map{meta, crams, crais, intervals -> [meta, crams, crais, intervals]},
        fasta,
        fai,
        dict,
        normal_vcf,
    )

    CLAIRS_PAIRED_FIRST.out.vcf_normal.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }.set{ clairs_vcf_normal_branch }

    clairs_vcf_normal_branch.intervals
        .map{ meta, vcf_normal_germline, vcf_normal_pileup ->
            [groupKey(meta, meta.num_intervals), vcf_normal_germline, vcf_normal_pileup]
        }.groupTuple().set{clairs_vcf_normal_branch_grouped}

    // Merge normal germline VCF
    //Only when using intervals

    clairs_vcf_normal_branch_grouped.map{meta, vcf_normal_germline, vcf_normal_pileup -> 
        new_meta = [
                        id:meta.normal_id + '_germline',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
        [new_meta, vcf_normal_germline]
    }.set{ch_clairs_vcf_normal_germline}

    clairs_vcf_normal_branch_grouped.map{meta, vcf_normal_germline, vcf_normal_pileup -> 
        new_meta = [
                        id:meta.normal_id + '_pileup',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
        [new_meta, vcf_normal_pileup]
    }.set{ch_clairs_vcf_normal_pileup}


    MERGE_VCFS_NORMAL_GERMLINE_CLAIRS(ch_clairs_vcf_normal_germline, dict)

    clairs_vcf_normal_germline = Channel.empty().mix(
        MERGE_VCFS_NORMAL_GERMLINE_CLAIRS.out.vcf,
        clairs_vcf_normal_branch.no_intervals.map{ meta, vcf_normal_germline, vcf_normal_pileup ->
            new_meta = [
                        id:meta.normal_id + '_germline',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
            [new_meta, vcf_normal_germline]
        }
    )

    MERGE_VCFS_NORMAL_PILEUP_CLAIRS(ch_clairs_vcf_normal_pileup, dict)

    clairs_vcf_normal_pileup = Channel.empty().mix(
        MERGE_VCFS_NORMAL_PILEUP_CLAIRS.out.vcf,
        clairs_vcf_normal_branch.no_intervals.map{ meta, vcf_normal_germline, vcf_normal_pileup ->
            new_meta = [
                        id:meta.normal_id + '_pileup',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
            [new_meta, vcf_normal_pileup]
        }
    )
    normal_vcf_for_rest = clairs_vcf_normal_germline.map{meta, g -> [g]}.first() ? clairs_vcf_normal_germline.map{meta, g -> [g]}.first() : normal_vcf
    CLAIRS_PAIRED_REST(
        input_branch.rest,
        fasta,
        fai,
        dict,
        normal_vcf_for_rest,
    )

    clairs_paired_vcfs = Channel.empty().mix(
        CLAIRS_PAIRED_FIRST.out.vcfs,
        CLAIRS_PAIRED_REST.out.vcfs
    )
    clairs_paired_tbi = Channel.empty().mix(
        CLAIRS_PAIRED_FIRST.out.tbi,
        CLAIRS_PAIRED_REST.out.tbi
    )


    // Figure out if using intervals or no_intervals
    clairs_paired_vcfs.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_vcfs_branch }

    clairs_paired_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_tbi_branch }
    
    clairs_vcfs_branch.intervals
        .map{ meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup ->
            [groupKey(meta, meta.num_intervals),  vcf_paired, vcf_tumor_germline, vcf_tumor_pileup]
        }.groupTuple().set{clairs_vcfs_branch_grouped}


    clairs_vcfs_branch_grouped.map{meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup -> 
        new_meta = [
                        id:meta.tumor_id + '_germline',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
        [new_meta, vcf_tumor_germline]
    }.set{ch_clairs_vcf_tumor_germline}

    clairs_vcfs_branch_grouped.map{meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup -> 
        new_meta = [
                        id:meta.tumor_id + '_pileup',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
        [new_meta, vcf_tumor_pileup]
    }.set{ch_clairs_vcf_tumor_pileup}
    

    clairs_vcfs_branch.intervals
        .map{ meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup ->
            [meta,  vcf_paired]
        }.set{ch_clairs_vcf_paired_for_fixheader}

    // Merge tumor germline VCF
    // Only when using intervals
    MERGE_VCFS_TUMOR_GERMLINE_CLAIRS(ch_clairs_vcf_tumor_germline, dict)

    clairs_vcf_tumor_germline = Channel.empty().mix(
        MERGE_VCFS_TUMOR_GERMLINE_CLAIRS.out.vcf,
        clairs_vcfs_branch.no_intervals.map{ meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup ->
            new_meta = [
                        id:meta.tumor_id + '_germline',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
            [new_meta, vcf_tumor_germline]
        }
    )

    // Merge tumor pilep VCF
    //Only when using intervals
    MERGE_VCFS_TUMOR_PILEUP_CLAIRS(ch_clairs_vcf_tumor_pileup, dict)

    clairs_vcf_tumor_pileup = Channel.empty().mix(
        MERGE_VCFS_TUMOR_PILEUP_CLAIRS.out.vcf,
        clairs_vcfs_branch.no_intervals.map{ meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup ->
            new_meta = [
                        id:meta.tumor_id + '_pileup',
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id
                    ]
            [new_meta, vcf_tumor_pileup]
        }
    )

    // Merge somatic VCF
    // First fix headers
    FIX_VCFHEADER_CLAIRS(
        ch_clairs_vcf_paired_for_fixheader,
        vcf_header
    )
    FIX_VCFHEADER_CLAIRS.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{ clairs_fixed_vcf_paired_branch }

    //Only when using intervals
    MERGE_VCFS_PAIRED_CLAIRS(
        clairs_fixed_vcf_paired_branch.intervals.map{meta, vcf_paired ->
            new_meta = [
                        id:meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id,
                        variantcaller: "clairs",
                    ]
            [groupKey(new_meta, meta.num_intervals), vcf_paired]}.groupTuple(),
        dict
    )

    clairs_vcf = Channel.empty().mix(
        MERGE_VCFS_PAIRED_CLAIRS.out.vcf,
        clairs_vcfs_branch.no_intervals.map{meta, vcf_paired, vcf_tumor_germline, vcf_tumor_pileup ->
            new_meta = [
                        id:meta.tumor_id + "_vs_" + meta.normal_id,
                        normal_id:meta.normal_id,
                        num_intervals:meta.num_intervals,
                        patient:meta.patient,
                        sex:meta.sex,
                        tumor_id:meta.tumor_id,
                        variantcaller: "clairs",
                    ]
            [new_meta, vcf_paired]
            }
        )

    clairs_tbi = Channel.empty().mix(
        MERGE_VCFS_PAIRED_CLAIRS.out.tbi,
        clairs_tbi_branch.no_intervals)

    ch_versions = ch_versions.mix(MERGE_VCFS_PAIRED_CLAIRS.out.versions)
    ch_versions = ch_versions.mix(CLAIRS_PAIRED_FIRST.out.versions)

    emit:
    clairs_vcf                  = clairs_vcf
    clairs_tbi                  = clairs_tbi
    clairs_vcf_tumor_germline   = clairs_vcf_tumor_germline
    clairs_tbi_tumor_germline   = MERGE_VCFS_TUMOR_GERMLINE_CLAIRS.out.tbi
    clairs_vcf_tumor_pileup     = clairs_vcf_tumor_pileup
    clairs_tbi_tumor_pileup     = MERGE_VCFS_TUMOR_PILEUP_CLAIRS.out.tbi 
    clairs_vcf_normal_germline  = clairs_vcf_normal_germline// channel: [ val(meta), [ vcf ] ]
    clairs_tbi_normal_germline  = MERGE_VCFS_NORMAL_GERMLINE_CLAIRS.out.tbi
    clairs_vcf_normal_pileup = clairs_vcf_normal_pileup// channel: [ val(meta), [ vcf ] ]
    clairs_tbi_normal_pileup  = MERGE_VCFS_NORMAL_PILEUP_CLAIRS.out.tbi

    versions                    = ch_versions                                   // channel: [ versions.yml ]
}
