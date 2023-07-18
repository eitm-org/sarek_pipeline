process MODKIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "$params.outdir/modkit", mode: 'symlink'
    container "ghcr.io/eitm-org/modkit"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path('*.bed')

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    bed = "${prefix}.bed"

    """
    echo "MODKIT $bam $bed"
    modkit pileup $bam $bed
    """
}
