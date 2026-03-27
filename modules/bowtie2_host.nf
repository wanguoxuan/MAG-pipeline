process BOWTIE2_HOST_REMOVE {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/host_removal", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(read1), path(read2)
    path host_index_dir

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_clean_R1.fq.gz"), path("${sample_id}_clean_R2.fq.gz"), emit: reads
    path "${sample_id}_host_removal_summary.txt", emit: summary

    script:
    """
    # Find the index basename from the .1.bt2 file
    bt2_file=\$(ls ${host_index_dir}/*.1.bt2 ${host_index_dir}/*.1.bt2l 2>/dev/null | head -1)
    index_base=\$(basename "\${bt2_file}" | sed 's/\\.1\\.bt2l\\?//')
    index_path="${host_index_dir}/\${index_base}"

    bowtie2 \
        -x \${index_path} \
        -1 ${read1} \
        -2 ${read2} \
        --un-conc-gz ${sample_id}_clean_R%.fq.gz \
        -p ${task.cpus} \
        --very-sensitive \
        -S /dev/null \
        2> ${sample_id}_host_removal_summary.txt
    """
}
