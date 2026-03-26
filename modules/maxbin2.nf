process MAXBIN2 {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/maxbin2", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs), path(depth_txt)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_maxbin2_bins"), emit: bins_dir
    path "${sample_id}_maxbin2.log",                                       emit: log

    script:
    """
    # Convert JGI depth file to MaxBin2 abundance format (contig_name <tab> mean_depth)
    awk 'NR > 1 { print \$1 "\\t" \$3 }' ${depth_txt} > ${sample_id}_abund.txt

    mkdir -p ${sample_id}_maxbin2_bins

    run_MaxBin.pl \
        -contig ${contigs} \
        -abund ${sample_id}_abund.txt \
        -out ${sample_id}_maxbin2_bins/${sample_id}_bin \
        -min_contig_length ${params.maxbin2_min_contig} \
        -thread ${task.cpus} \
        2>&1 | tee ${sample_id}_maxbin2.log

    # Rename .fasta extension to .fa for consistency
    for f in ${sample_id}_maxbin2_bins/*.fasta; do
        [ -e "\$f" ] && mv "\$f" "\${f%.fasta}.fa"
    done
    """
}
