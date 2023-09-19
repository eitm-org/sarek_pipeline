
process MODKIT {

    tag "$meta.id"
    label 'process_high'
    container "ghcr.io/eitm-org/modkit:1.0"

    input:

    tuple val(meta), path(bam), path(bai)
    path fasta

    output:

    path('*.bed')
    path(summary)
    path(modkit_log)
    path("versions.yml"), emit: versions

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def bed = "${prefix}.bed"
    summary = "${prefix}.summary"
    modkit_log = "${prefix}.log"

    """
    echo "MODKIT $bam $bed"

    date >> $modkit_log
    echo "modkit summary" >> $modkit_log
    modkit summary $bam --log-filepath $modkit_log > $summary 

    date >> $modkit_log
    echo "modkit pileup" >> $modkit_log
    echo "args: $args" >> $modkit_log
    modkit pileup $args $bam $bed --log-filepath $modkit_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(echo \$(modkit --version) | sed 's/mod_kit //')
    END_VERSIONS
    """
}
