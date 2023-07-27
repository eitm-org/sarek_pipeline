
process MODKIT {

    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/eitm-org/modkit"

    publishDir "$params.outdir/modkit", mode: params.publish_dir_mode
    // TODO: do we want to move publishDir info to conf/modules/modules.config

    input:

    tuple val(meta), path(bam), path(bai)
    path fasta

    output:

    path('*.bed')
    path(summary)
    path(modkit_log)
    path("versions.yml"), emit: versions

    script:

    prefix = task.ext.prefix ?: "${meta.id}"
    bed = "${prefix}.bed"
    summary = "${prefix}.summary"
    modkit_log = "${prefix}.log"

    // modkit command line options:
    //  https://nanoporetech.github.io/modkit/

    modkit_args = " --ref $fasta"
    modkit_args += " -t $task.cpus"
    modkit_args += " --preset traditional" // "--cpg --ignore h --combine-strands"

    // TODO: investigate --interval-size / --region / --include-bed 
    // for parallelization performance

    """
    echo "MODKIT $bam $bed"

    date >> $modkit_log
    echo "modkit summary" >> $modkit_log
    modkit summary $bam --log-filepath $modkit_log > $summary 

    date >> $modkit_log
    echo "modkit pileup" >> $modkit_log
    echo "modkit_args: $modkit_args" >> $modkit_log
    modkit pileup $modkit_args $bam $bed --log-filepath $modkit_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(echo \$(modkit --version) | sed 's/mod_kit //')
    END_VERSIONS
    """
}
