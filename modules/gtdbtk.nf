process GTDBTK {
    tag "${sample_id}"
    label 'process_gtdbtk'
    publishDir "${params.output_dir}/${batch}/${sample_id}/gtdbtk", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(bins_dir)
    val gtdbtk_db

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_gtdbtk"),                    emit: classify
    path "${sample_id}_gtdbtk/gtdbtk.bac120.summary.tsv", optional: true,             emit: bac_summary
    path "${sample_id}_gtdbtk/gtdbtk.ar53.summary.tsv",   optional: true,             emit: arc_summary

    script:
    """
    # Skip if no bins were produced
    n_bins=\$(ls ${bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$n_bins" -eq 0 ]; then
        echo "No bins found for ${sample_id}, skipping GTDB-Tk" >&2
        mkdir -p ${sample_id}_gtdbtk
        exit 0
    fi

    export GTDBTK_DATA_PATH=${gtdbtk_db}

    gtdbtk classify_wf \
        --genome_dir ${bins_dir} \
        --out_dir ${sample_id}_gtdbtk \
        --cpus ${task.cpus} \
        --extension fa \
        ${params.gtdbtk_skani_sketch_dir ? "--skani_sketch_dir ${params.gtdbtk_skani_sketch_dir}" : "--skip_ani_screen"}
    """
}

process MERGE_GTDBTK_SUMMARIES {
    label 'process_low'
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    path bac_summaries
    path arc_summaries

    output:
    path "all_gtdbtk_bac120_summary.tsv", optional: true, emit: bac_merged
    path "all_gtdbtk_ar53_summary.tsv",   optional: true, emit: arc_merged

    script:
    """
    # Merge bacterial summaries
    bac_files=(\$(ls *bac120.summary.tsv 2>/dev/null || true))
    if [ \${#bac_files[@]} -gt 0 ]; then
        head -1 "\${bac_files[0]}" > all_gtdbtk_bac120_summary.tsv
        for f in "\${bac_files[@]}"; do
            tail -n +2 "\$f" >> all_gtdbtk_bac120_summary.tsv
        done
    fi

    # Merge archaeal summaries
    arc_files=(\$(ls *ar53.summary.tsv 2>/dev/null || true))
    if [ \${#arc_files[@]} -gt 0 ]; then
        head -1 "\${arc_files[0]}" > all_gtdbtk_ar53_summary.tsv
        for f in "\${arc_files[@]}"; do
            tail -n +2 "\$f" >> all_gtdbtk_ar53_summary.tsv
        done
    fi
    """
}
