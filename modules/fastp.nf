process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.output_dir}/${batch}/${sample_id}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(read1), path(read2)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_R1_trimmed.fq.gz"), path("${sample_id}_R2_trimmed.fq.gz"), emit: reads
    path "${sample_id}_fastp.json", emit: json
    path "${sample_id}_fastp.html", emit: html

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}_R1_trimmed.fq.gz \
        --out2 ${sample_id}_R2_trimmed.fq.gz \
        --json ${sample_id}_fastp.json \
        --html ${sample_id}_fastp.html \
        --qualified_quality_phred ${params.fastp_min_qual} \
        --length_required ${params.fastp_min_len} \
        --detect_adapter_for_pe \
        --thread ${task.cpus}
    """
}
