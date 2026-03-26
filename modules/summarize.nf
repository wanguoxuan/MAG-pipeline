process SUMMARIZE_BINS {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.output_dir}/${batch}/${sample_id}/metabat2", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(bins_dir)

    output:
    path "${sample_id}_bin_summary.tsv", emit: summary

    script:
    """
    echo -e "sample_id\tbin_id\tcontigs\ttotal_length_bp" > ${sample_id}_bin_summary.tsv

    shopt -s nullglob
    for bin in ${bins_dir}/*.fa; do
        bin_id=\$(basename \$bin .fa)
        n_contigs=\$(grep -c "^>" \$bin)
        total_len=\$(grep -v "^>" \$bin | tr -d '\\n' | wc -c)
        echo -e "${sample_id}\\t\${bin_id}\\t\${n_contigs}\\t\${total_len}" >> ${sample_id}_bin_summary.tsv
    done

    # Report if no bins were produced
    if [ \$(wc -l < ${sample_id}_bin_summary.tsv) -eq 1 ]; then
        echo -e "${sample_id}\\tno_bins\\t0\\t0" >> ${sample_id}_bin_summary.tsv
    fi
    """
}

process MERGE_BIN_SUMMARIES {
    label 'process_low'
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    path summaries

    output:
    path "all_bins_summary.tsv", emit: merged

    script:
    """
    # Write header once, then append data rows from all samples
    head -1 \$(ls *_bin_summary.tsv | head -1) > all_bins_summary.tsv
    for f in *_bin_summary.tsv; do
        tail -n +2 \$f >> all_bins_summary.tsv
    done
    """
}
