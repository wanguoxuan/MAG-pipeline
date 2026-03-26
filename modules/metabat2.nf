process JGI_DEPTH {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.output_dir}/${batch}/${sample_id}/metabat2", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs), path(bam), path(bai)

    output:
    tuple val(sample_id), val(batch), path(contigs), path("${sample_id}_depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths \
        --outputDepth ${sample_id}_depth.txt \
        ${bam}
    """
}

process METABAT2 {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/metabat2", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs), path(depth)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_metabat2_bins"), emit: bins
    path "${sample_id}_metabat2.log",                                      emit: log

    script:
    """
    mkdir -p ${sample_id}_metabat2_bins
    metabat2 \
        -i ${contigs} \
        -a ${depth} \
        -o ${sample_id}_metabat2_bins/${sample_id}_bin \
        --minContig ${params.metabat2_min_contig} \
        -t ${task.cpus} \
        -v \
        2> ${sample_id}_metabat2.log
    """
}
