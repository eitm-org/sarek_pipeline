process MODKIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "$params.outdir/modkit", mode: 'symlink'
    container "ghcr.io/eitm-org/modkit"
    //containerOptions "-v $launchDir/$params.outdir/modkit:/modkit_output"
    time '1h'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    //tuple path(bam), path('*.bed')
    path('*.bed')

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    //bed = "/modkit_output/${bam.toString().replace(/.bam/, /.bed/)}"
    bed = "${prefix}.bed"

    """
    echo "MODKIT $bam $bed"
    touch $bed
    modkit pileup $bam $bed
    """
}
