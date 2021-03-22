// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_BAM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::samtools=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.11--h6270b1f_0"
    } else {
        container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${name}.bam"), emit: bam
        val meta,                            emit: tsv

    script:
    name = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    samtools merge --threads ${task.cpus} ${name}.bam ${bam}
    """
}
