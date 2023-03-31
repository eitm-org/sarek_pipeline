process GATK4_FIXVCFHEADER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path  vcf_header

    output:
    tuple val(meta), path('*_fixedheader.vcf.gz'), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"


    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MergeVcfs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" FixVcfHeader \\
        -I $vcf \\
        --OUTPUT tmp.vcf.gz\\
        --HEADER $vcf_header \\
        $args
    
    gatk --java-options "-Xmx${avail_mem}g" RenameSampleInVcf \\
        --INPUT tmp.vcf.gz\\
        --OUTPUT ${vcf.baseName}_fixedheader.vcf.gz \\
        --NEW_SAMPLE_NAME ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${vcf.baseName}_fixedheader.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
