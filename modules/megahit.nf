process MEGAHIT {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/megahit", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(read1), path(read2)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_contigs.fa"), emit: contigs
    path "${sample_id}_megahit.log",                                    emit: log

    script:
    """
    megahit \
        -1 ${read1} \
        -2 ${read2} \
        --min-contig-len ${params.min_contig_len} \
        --num-cpu-threads ${task.cpus} \
        -o megahit_out \
        --out-prefix ${sample_id} \
        2> ${sample_id}_megahit.log

    cp megahit_out/${sample_id}.contigs.fa ${sample_id}_contigs.fa
    """
}
