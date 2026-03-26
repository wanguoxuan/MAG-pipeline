#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * MAG Binning Pipeline
 * fastp → Host Removal (optional) → MEGAHIT → Bowtie2 → MetaBat2 + MaxBin2 + CONCOCT → DAS Tool → CheckM2 → GTDB-Tk
 */

// ========================
// Parameter Definitions
// ========================

params.sample_sheet  = null   // CSV with columns: sample_id,batch,read1,read2
params.bins_sheet    = null   // CSV with columns: sample_id,batch,bins_dir (for GTDBTK_ONLY entry)
params.output_dir    = "./results"

// QC parameters
params.fastp_min_qual  = 20
params.fastp_min_len   = 50
params.fastp_cpus      = 8

// Host removal (optional — set to null to skip)
params.host_index      = null   // Path to directory containing Bowtie2 index for host genome
params.bowtie2_cpus    = 24

// Assembly parameters
params.min_contig_len      = 500
params.megahit_cpus        = 24

// Binning parameters
params.metabat2_min_contig  = 1500
params.metabat2_cpus        = 16
params.maxbin2_min_contig   = 1500
params.maxbin2_cpus         = 16
params.concoct_chunk_size   = 10000
params.concoct_min_contig   = 1000
params.concoct_cpus         = 16

// Bin refinement
params.dastool_score_threshold = 0.5
params.dastool_cpus            = 16

// Quality check
params.checkm2_db    = null   // Path to CheckM2 database (optional; uses CHECKM2DB env var if null)
params.checkm2_cpus  = 16

// Taxonomy
params.gtdbtk_db     = null
params.gtdbtk_cpus   = 24

// Skip options
params.skip_summary  = false
params.skip_gtdbtk   = false

// Help
params.help = false

def helpMessage() {
    log.info"""
    =========================================
    MAG Binning Pipeline
    =========================================
    Usage:
      nextflow run main.nf --sample_sheet samples.csv --output_dir results

    Required:
      --sample_sheet          CSV file with columns: sample_id,batch,read1,read2
      --output_dir            Output directory [${params.output_dir}]

    QC & Trimming (fastp):
      --fastp_min_qual        Minimum base quality [${params.fastp_min_qual}]
      --fastp_min_len         Minimum read length [${params.fastp_min_len}]
      --fastp_cpus            CPUs for fastp [${params.fastp_cpus}]

    Host Removal (optional):
      --host_index            Path to directory containing Bowtie2 host genome index
                              (skip host removal if not provided)
      --bowtie2_cpus          CPUs for Bowtie2 [${params.bowtie2_cpus}]

    Assembly (MEGAHIT):
      --min_contig_len        Minimum contig length to keep [${params.min_contig_len}]
      --megahit_cpus          CPUs for MEGAHIT [${params.megahit_cpus}]

    Binning:
      --metabat2_min_contig   MetaBat2 minimum contig length [${params.metabat2_min_contig}]
      --metabat2_cpus         CPUs for MetaBat2 [${params.metabat2_cpus}]
      --maxbin2_min_contig    MaxBin2 minimum contig length [${params.maxbin2_min_contig}]
      --maxbin2_cpus          CPUs for MaxBin2 [${params.maxbin2_cpus}]
      --concoct_chunk_size    CONCOCT contig chunk size [${params.concoct_chunk_size}]
      --concoct_min_contig    CONCOCT minimum contig length [${params.concoct_min_contig}]
      --concoct_cpus          CPUs for CONCOCT [${params.concoct_cpus}]

    Bin Refinement (DAS Tool):
      --dastool_score_threshold  DAS Tool score threshold [${params.dastool_score_threshold}]
      --dastool_cpus             CPUs for DAS Tool [${params.dastool_cpus}]

    Quality Check (CheckM2):
      --checkm2_db            Path to CheckM2 database (uses CHECKM2DB env var if not set)
      --checkm2_cpus          CPUs for CheckM2 [${params.checkm2_cpus}]

    Taxonomy (GTDB-Tk):
      --gtdbtk_db             Path to GTDB-Tk reference data (required unless --skip_gtdbtk)
      --gtdbtk_cpus           CPUs for GTDB-Tk [${params.gtdbtk_cpus}]
      --skip_gtdbtk           Skip GTDB-Tk taxonomic classification

    Other:
      --skip_summary          Skip summary generation
      --help                  Show this message
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Entry-specific input validation is done inside each workflow block

if (!params.skip_gtdbtk && !params.gtdbtk_db) {
    log.error "ERROR: --gtdbtk_db is required unless --skip_gtdbtk is set"
    exit 1
}

// ========================
// Include Modules
// ========================

include { FASTP }                  from './modules/fastp'
include { BOWTIE2_HOST_REMOVE }    from './modules/bowtie2_host'
include { MEGAHIT }                from './modules/megahit'
include { BOWTIE2_MAP }            from './modules/bowtie2_mapping'
include { JGI_DEPTH }              from './modules/metabat2'
include { METABAT2 }               from './modules/metabat2'
include { MAXBIN2 }                from './modules/maxbin2'
include { CONCOCT }                from './modules/concoct'
include { DASTOOL }                from './modules/dastool'
include { CHECKM2 }                from './modules/checkm2'
include { MERGE_CHECKM2_SUMMARIES } from './modules/checkm2'
include { SUMMARIZE_BINS }         from './modules/summarize'
include { MERGE_BIN_SUMMARIES }    from './modules/summarize'
include { GTDBTK }                 from './modules/gtdbtk'
include { MERGE_GTDBTK_SUMMARIES } from './modules/gtdbtk'

// ========================
// Main Workflow
// ========================

// ========================
// GTDB-Tk only entry point
// ========================

workflow GTDBTK_ONLY {

    if (!params.bins_sheet) {
        log.error "ERROR: --bins_sheet is required for GTDBTK_ONLY entry (CSV with columns: sample_id,batch,bins_dir)"
        exit 1
    }

    Channel
        .fromPath(params.bins_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.batch, file(row.bins_dir)) }
        .set { bins_ch }

    GTDBTK(bins_ch, params.gtdbtk_db)

    MERGE_GTDBTK_SUMMARIES(
        GTDBTK.out.bac_summary.collect(),
        GTDBTK.out.arc_summary.collect()
    )
}

// ========================
// Full pipeline entry point
// ========================

workflow {

    if (!params.sample_sheet) {
        log.error "ERROR: --sample_sheet is required"
        helpMessage()
        exit 1
    }

    // Read sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                row.batch,
                file(row.read1),
                file(row.read2)
            )
        }
        .set { samples_ch }

    // Step 1: QC & trimming with fastp
    FASTP(samples_ch)

    // Step 2: Optional host read removal
    if (params.host_index) {
        host_index_ch = Channel.fromPath(params.host_index, checkIfExists: true)
        BOWTIE2_HOST_REMOVE(FASTP.out.reads, host_index_ch)
        clean_reads_ch = BOWTIE2_HOST_REMOVE.out.reads
    } else {
        clean_reads_ch = FASTP.out.reads
    }

    // Step 3: Assemble reads into contigs with MEGAHIT
    MEGAHIT(clean_reads_ch)

    // Step 4: Map reads back to assembly for coverage depth
    MEGAHIT.out.contigs
        .join(clean_reads_ch.map { sid, batch, r1, r2 -> tuple(sid, r1, r2) })
        .map { sid, batch, contigs, r1, r2 -> tuple(sid, batch, contigs, r1, r2) }
        .set { mapping_input_ch }

    BOWTIE2_MAP(mapping_input_ch)

    // Step 5: Calculate per-contig depth from BAM
    JGI_DEPTH(BOWTIE2_MAP.out.bam)
    // JGI_DEPTH.out.depth: (sample_id, batch, contigs, depth_txt)

    // Step 6: Parallel binning with three binners
    METABAT2(JGI_DEPTH.out.depth)
    MAXBIN2(JGI_DEPTH.out.depth)

    // CONCOCT needs contigs + BAM; extract from existing channels
    JGI_DEPTH.out.depth
        .map { sid, batch, contigs, _depth -> tuple(sid, batch, contigs) }
        .join(BOWTIE2_MAP.out.bam.map { sid, _batch, _ctg, bam, bai -> tuple(sid, bam, bai) })
        .map { sid, batch, contigs, bam, bai -> tuple(sid, batch, contigs, bam, bai) }
        .set { concoct_input_ch }

    CONCOCT(concoct_input_ch)

    // Step 7: Bin refinement with DAS Tool
    // Join binner outputs on sample_id and combine with contigs
    METABAT2.out.bins
        .join(MAXBIN2.out.bins_dir)
        .map { sid, b, mb2, _b2, mx2 -> tuple(sid, b, mb2, mx2) }
        .join(CONCOCT.out.bins_dir)
        .map { sid, b, mb2, mx2, _b3, co -> tuple(sid, b, mb2, mx2, co) }
        .join(JGI_DEPTH.out.depth.map { sid, b, contigs, _d -> tuple(sid, b, contigs) })
        .map { sid, b, mb2, mx2, co, _b4, ctg -> tuple(sid, b, ctg, mb2, mx2, co) }
        .set { dastool_input_ch }

    DASTOOL(dastool_input_ch)

    // Step 8: Quality check with CheckM2
    CHECKM2(DASTOOL.out.bins)

    if (!params.skip_summary) {
        SUMMARIZE_BINS(DASTOOL.out.bins)

        SUMMARIZE_BINS.out.summary
            .collect()
            .set { all_bin_summaries_ch }

        MERGE_BIN_SUMMARIES(all_bin_summaries_ch)

        CHECKM2.out.quality_report
            .collect()
            .set { all_checkm2_ch }

        MERGE_CHECKM2_SUMMARIES(all_checkm2_ch)
    }

    // Step 9: Taxonomic classification with GTDB-Tk
    if (!params.skip_gtdbtk) {
        GTDBTK(DASTOOL.out.bins, params.gtdbtk_db)

        MERGE_GTDBTK_SUMMARIES(
            GTDBTK.out.bac_summary.collect(),
            GTDBTK.out.arc_summary.collect()
        )
    }
}

workflow.onComplete {
    log.info """
    =========================================
    MAG Binning Pipeline completed!
    =========================================
    Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:  ${workflow.duration}
    Output:    ${params.output_dir}
    =========================================
    """.stripIndent()
}
