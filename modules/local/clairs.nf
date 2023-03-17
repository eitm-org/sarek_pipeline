process CLAIRS {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs:latest':
        'hkubal/clairs:latest' }"

    input:
    tuple val(meta), path(input_normal), path(input_normal_index), path(intervals)
    tuple val(meta), path(input_tumor), path(input_tumor_index), path(intervals)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    what
    def args = task.ext.args ?: ''
    def inputs_normal = input_normal.collect{ "--normal_bam_fn $it"}.join(" ")
    def inputs_tumor = input_normal.collect{ "--tumor_bam_fn $it"}.join(" ")
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region_command = intervals ? "--region $intervals" : ""

    """
    /opt/bin/run_clairs \\
        --normal_bam_fn ${input.normal_cram} \\
        --tumor_bam_fn ${input.tumor_cram} \\
        --ref_fn ${fasta} \\
        --threads ${task.cpus} \\
        --platform ont_r10 \\
        --output_dir . \\
        --output_prefix $prefix \\
        $region_command \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(echo \$(/opt/bin/run_clairs --version 2>&1) | sed 's/^.*(clairS) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch run_clairs.log.bak
    touch run_clairs.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
