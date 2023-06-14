//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../../modules/nf-core/samtools/merge/main'


process SORT_BAM {
    
    input:
    tuple val(meta), path(input_files, stageAs: "?/*")
 
    // TODO: remove stageAs?

    output:
    tuple val(meta), path(output_files), optional:true, emit: bam

    // TODO: remove optional, emit?

    script:

    // BAM files are staged in subdirs of the work directory:
    //   1/file_234.bam 
    //   2/file_123.bam
    //   ...
    // We want to sort based on the filename, so we remove
    // the relative path to the bam file for the sort criterion.

    // this works, but danger of array out of bounds
    //output_files = input_files.sort {it.name.split('/')[1]}

    output_files = input_files.sort {it.name.replaceAll('.*/','')}

    println "***** input_files: " + input_files
    println "***** output_files: " + output_files

    """
    echo '***** sorting'
    """
}


workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
        bam // channel: [mandatory] meta, bam

    main:
    ch_versions = Channel.empty()

    // Figuring out if there is one or more bam(s) from the same sample
    bam.branch{
        //Here there actually is a list, so size() works
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge}

    // TODO dk: sort before merge
    //

    // TODO: remove
    bam_to_merge.single.view { "***** bam_to_merge.single: $it" }
    bam_to_merge.multiple.view { "***** bam_to_merge.multiple: $it" }

    bam_to_merge_sorted = SORT_BAM(bam_to_merge.multiple)

    // TODO: remove
    bam_to_merge_sorted.view { "***** bam_to_merge_sorted: $it" }
    SORT_BAM.out.bam.view { "***** SORT_BAM.out.bam: $it" }

    MERGE_BAM(bam_to_merge_sorted, [], [])
    INDEX_MERGE_BAM(bam_to_merge.single.mix(MERGE_BAM.out.bam))

    bam_bai = bam_to_merge.single.map{meta, bam -> [meta, bam[0]]}
        .mix(MERGE_BAM.out.bam)
        .join(INDEX_MERGE_BAM.out.bai)

    //dk
    bam_bai.view { "***** bam_bai: $it" }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_MERGE_BAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())

    emit:
        bam_bai
        versions = ch_versions
}
