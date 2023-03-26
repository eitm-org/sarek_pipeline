process CLAIRS {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs:latest':
        'clairs_input_normal_vcf:latest' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    path fasta
    path fai
    path dict
    path(normal_vcf)

    output:
    tuple   val(meta), 
            path("*.clairs.*.vcf.gz"),
            path("*_tumor_germline_*.vcf.gz"),
            path("*_tumor_pileup_*.vcf.gz"),             emit: vcfs
    tuple val(meta), path("*_normal_germline_*.vcf.gz"), path("*_normal_pileup_*.vcf.gz"), emit: vcf_normal, optional: true
    tuple val(meta), path("*.clairs.*.vcf.gz.tbi"),     emit: tbi
    path "versions.yml",                                emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    
    def args = task.ext.args ?: ''
    def inputs = "--normal_bam_fn ${input[0]} --tumor_bam_fn ${input[1]}"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = intervals ? "${intervals.toString()}" : ""
    def bed_command = intervals ? "-b ${intervals}" : ""
    def normal_vcf_fn = normal_vcf ? "--normal_vcf_fn ${normal_vcf}" : "" 
    def mv_normal_vcf_germline_command = normal_vcf ? "" : "mv tmp/clair3_output/clair3_normal_output/merge_output.vcf.gz ${meta.normal_id}_normal_germline_${suffix}.vcf.gz"
    def mv_normal_vcf_pileup_command = normal_vcf ? "" : "mv tmp/clair3_output/clair3_normal_output/pileup.vcf.gz ${meta.normal_id}_normal_pileup_${suffix}.vcf.gz"

    """
    /opt/bin/run_clairs \\
        $inputs \\
        $normal_vcf_fn \\
        --ref_fn ${fasta} \\
        --threads ${task.cpus} \\
        --output_dir . \\
        --output_prefix $prefix \\
        $bed_command \\
        $args


    $mv_normal_vcf_germline_command
    $mv_normal_vcf_pileup_command

    mv tmp/clair3_output/clair3_tumor_output/merge_output.vcf.gz ${meta.tumor_id}_tumor_germline_${suffix}.vcf.gz
    mv tmp/clair3_output/clair3_tumor_output/pileup.vcf.gz ${meta.tumor_id}_tumor_pileup_${suffix}.vcf.gz

    
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
    touch ${meta.normal_id}_normal_germline_${suffix}.vcf.gz
    touch ${meta.normal_id}_normal_pileup_${suffix}.vcf.gz
    touch ${meta.tumor_id}_tumor_germline_${suffix}.vcf.gz
    touch ${meta.tumor_id}_tumor_pileup_${suffix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
