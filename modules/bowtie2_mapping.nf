process BOWTIE2_MAP {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/mapping", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs), path(read1), path(read2)

    output:
    tuple val(sample_id), val(batch), path(contigs), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam
    path "${sample_id}_mapping_summary.txt", emit: summary

    script:
    """
    # Build Bowtie2 index from assembly contigs
    bowtie2-build \
        --threads ${task.cpus} \
        ${contigs} \
        contigs_index

    # Map reads back to contigs and sort BAM
    bowtie2 \
        -p ${task.cpus} \
        -x contigs_index \
        -1 ${read1} \
        -2 ${read2} \
        --very-sensitive \
        2> ${sample_id}_mapping_summary.txt \
    | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam

    samtools index ${sample_id}_sorted.bam
    """
}
