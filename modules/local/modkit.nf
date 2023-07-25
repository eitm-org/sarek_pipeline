process MODKIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "$params.outdir/modkit", mode: 'symlink'
    container "ghcr.io/eitm-org/modkit"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path('*.bed')
    path(summary)
    path(modkit_log)

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    bed = "${prefix}.bed"
    summary = "${prefix}.summary"
    modkit_log = "${prefix}.log"

    """
    echo "MODKIT $bam $bed"

    date >> $modkit_log
    echo "modkit summary" >> $modkit_log
    modkit summary $bam --log-filepath $modkit_log > $summary 

    date >> $modkit_log
    echo "modkit pileup" >> $modkit_log
    modkit pileup $bam $bed --log-filepath $modkit_log
    """
}
