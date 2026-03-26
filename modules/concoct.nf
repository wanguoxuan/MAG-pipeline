process CONCOCT {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.output_dir}/${batch}/${sample_id}/concoct", mode: 'copy'

    input:
    tuple val(sample_id), val(batch), path(contigs), path(bam), path(bai)

    output:
    tuple val(sample_id), val(batch), path("${sample_id}_concoct_bins"), emit: bins_dir
    path "${sample_id}_concoct_output",                                   emit: output_dir

    script:
    """
    # Cut contigs into chunks for coverage calculation
    cut_up_fasta.py ${contigs} \
        -c ${params.concoct_chunk_size} \
        -o 0 \
        --merge_last \
        -b ${sample_id}_contigs_10k.bed \
        > ${sample_id}_contigs_10k.fa

    # Generate coverage table from BAM
    concoct_coverage_table.py \
        ${sample_id}_contigs_10k.bed \
        ${bam} \
        > ${sample_id}_coverage_table.tsv

    # Run CONCOCT clustering
    concoct \
        --composition_file ${sample_id}_contigs_10k.fa \
        --coverage_file ${sample_id}_coverage_table.tsv \
        -b ${sample_id}_concoct_output/ \
        --threads ${task.cpus} \
        -l ${params.concoct_min_contig}

    # Merge subcontig cluster assignments back to original contigs
    merge_cutup_clustering.py \
        ${sample_id}_concoct_output/clustering_gt${params.concoct_min_contig}.csv \
        > ${sample_id}_concoct_output/clustering_merged.csv

    # Extract bins as FASTA files
    mkdir -p ${sample_id}_concoct_bins
    extract_fasta_bins.py \
        ${contigs} \
        ${sample_id}_concoct_output/clustering_merged.csv \
        --output_path ${sample_id}_concoct_bins

    # Rename .fa extension for consistency (CONCOCT outputs .fa by default)
    for f in ${sample_id}_concoct_bins/*.fa; do
        true  # already .fa
    done
    """
}
