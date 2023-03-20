process CLAIRS {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs:latest':
        'hkubal/clairs:latest' }"

    input:
    tuple val(meta), path(input), path(inputl_index), path(intervals)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    def args = task.ext.args ?: ''
    def inputs = "--normal_bam_fn ${input[0]} --tumor_bam_fn ${input[1]}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def region_command = intervals ? "-c ${intervals.toString().split('_')[0]}" : ""

    """
    /opt/bin/run_clairs \\
        $inputs \\
        --ref_fn ${fasta} \\
        --threads ${task.cpus} \\
        --platform ont_r10 \\
        --output_dir . \\
        --output_prefix $prefix \\
        --disable_phasing \\
        $region_command \\
        $args
    
    mv tmp/* ..
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
