process CHECKM2 {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(bins_dir)

    output:
    tuple val(sample_id), val(batch), path("checkm2"),                           emit: output_dir
    path "checkm2/quality_report.tsv", optional: true,                            emit: quality_report

    script:
    """
    n_bins=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$n_bins" -eq 0 ]; then
        echo "No bins found for ${sample_id}, skipping CheckM2" >&2
        mkdir -p checkm2
        exit 0
    fi

    checkm2 predict \
        --input ${bins_dir} \
        --output-directory checkm2 \
        --extension fa \
        --threads ${task.cpus} \
        ${params.checkm2_db ? "--database_path ${params.checkm2_db}" : ""}
    """
}

process MERGE_CHECKM2_SUMMARIES {
    label 'process_low'
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    path quality_reports

    output:
    path "all_checkm2_quality.tsv", emit: merged

    script:
    """
    reports=(\$(ls *quality_report.tsv 2>/dev/null || true))
    if [ \${#reports[@]} -eq 0 ]; then
        echo "No CheckM2 reports found" >&2
        exit 0
    fi
    head -1 "\${reports[0]}" > all_checkm2_quality.tsv
    for f in "\${reports[@]}"; do
        tail -n +2 "\$f" >> all_checkm2_quality.tsv
    done
    """
}
