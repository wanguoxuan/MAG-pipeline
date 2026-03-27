process DASTOOL {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/dastool", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs),
          path(metabat2_bins_dir), path(maxbin2_bins_dir), path(concoct_bins_dir)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_dastool_bins"), emit: bins
    path "${sample_id}_dastool_DASTool_summary.tsv", optional: true,     emit: summary
    path "${sample_id}_dastool_DASTool_scores.tsv",  optional: true,     emit: scores

    script:
    """
    s2b_files=""
    labels=""

    # Build scaffold2bin for MetaBat2 bins
    n_metabat=\$(ls ${metabat2_bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$n_metabat" -gt 0 ]; then
        Fasta_to_Contig2Bin.sh -e fa -i ${metabat2_bins_dir} > metabat2_s2b.tsv
        s2b_files="metabat2_s2b.tsv"
        labels="metabat2"
    fi

    # Build scaffold2bin for MaxBin2 bins
    n_maxbin=\$(ls ${maxbin2_bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$n_maxbin" -gt 0 ]; then
        Fasta_to_Contig2Bin.sh -e fa -i ${maxbin2_bins_dir} > maxbin2_s2b.tsv
        s2b_files="\${s2b_files:+\${s2b_files},}maxbin2_s2b.tsv"
        labels="\${labels:+\${labels},}maxbin2"
    fi

    # Build scaffold2bin for CONCOCT bins
    n_concoct=\$(ls ${concoct_bins_dir}/*.fa 2>/dev/null | wc -l)
    if [ "\$n_concoct" -gt 0 ]; then
        Fasta_to_Contig2Bin.sh -e fa -i ${concoct_bins_dir} > concoct_s2b.tsv
        s2b_files="\${s2b_files:+\${s2b_files},}concoct_s2b.tsv"
        labels="\${labels:+\${labels},}concoct"
    fi

    # Exit if no bins from any binner
    if [ -z "\$s2b_files" ]; then
        echo "No bins from any binner for ${sample_id}, skipping DAS Tool" >&2
        mkdir -p ${sample_id}_dastool_bins
        exit 0
    fi

    # Run DAS Tool
    DAS_Tool \
        -i "\$s2b_files" \
        -l "\$labels" \
        -c ${contigs} \
        -o ${sample_id}_dastool \
        --score_threshold ${params.dastool_score_threshold} \
        --threads ${task.cpus} \
        --write_bins

    # Move refined bins to output directory
    if [ -d "${sample_id}_dastool_DASTool_bins" ]; then
        mv ${sample_id}_dastool_DASTool_bins ${sample_id}_dastool_bins
    else
        mkdir -p ${sample_id}_dastool_bins
    fi
    """
}
