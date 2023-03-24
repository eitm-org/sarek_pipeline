process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--haef29d1_2':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"
        // 'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path regions
    path targets
    path samples

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions ? "--regions-file ${regions}" : ""
    def targets_file = targets ? "--targets-file ${targets}" : ""
    def samples_file =  samples ? "--samples-file ${samples}" : ""
    def sample_command = tbi ? "-s SAMPLE" : ""
    """
    bcftools +fill-tags $vcf -o $vcf -Ov -- -t all

    bcftools stats \\
        --verbose \\
        $args \\
        $targets_file \\
        $samples_file \\
        $sample_command \\
        $regions_file \\
        $vcf > ${prefix}.bcftools_stats.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
