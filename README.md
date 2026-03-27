# MAG Binning Pipeline

A Nextflow pipeline for Metagenome-Assembled Genome (MAG) binning and taxonomic classification.

## Workflow

```
reads → FASTP → trimmed reads
trimmed reads → BOWTIE2_HOST_REMOVE (optional) → clean reads
clean reads → MEGAHIT → contigs
clean reads + contigs → BOWTIE2_MAP → sorted BAM
contigs + BAM → JGI_DEPTH → depth.txt
                                    ├── contigs + depth.txt → METABAT2 → metabat2_bins/
                                    ├── contigs + depth.txt → MAXBIN2  → maxbin2_bins/
                                    └── contigs + BAM       → CONCOCT  → concoct_bins/
metabat2_bins + maxbin2_bins + concoct_bins + contigs → DASTOOL → refined_bins/
refined_bins → CHECKM2 → quality_report.tsv
refined_bins → GTDBTK  → per-sample taxonomy TSVs
all taxonomy TSVs → MERGE_GTDBTK_SUMMARIES → all_gtdbtk_*.tsv
all quality TSVs  → MERGE_CHECKM2_SUMMARIES → all_checkm2_quality.tsv
```

## Tools

| Step | Tool | Version |
|------|------|---------|
| QC & Trimming | fastp | 1.0.1 |
| Host Removal | Bowtie2 + Samtools | 2.4.5 / 1.15.1 |
| Assembly | MEGAHIT | 1.2.9 |
| Read Mapping | Bowtie2 + Samtools | 2.4.5 / 1.15.1 |
| Depth calculation | jgi_summarize_bam_contig_depths | MetaBat2 2.15 |
| Binning | MetaBat2 | 2.15 |
| Binning | MaxBin2 | 2.2.7 |
| Binning | CONCOCT | 1.1.0 |
| Bin Refinement | DAS Tool | 1.1.6 |
| Quality Check | CheckM2 | 1.1.0 |
| Taxonomy | GTDB-Tk | 2.6.1 |

## Usage

### Required input

A sample sheet CSV with columns: `sample_id`, `batch`, `read1`, `read2`

```csv
sample_id,batch,read1,read2
S001,BATCH_1,/path/to/S001_R1.fq.gz,/path/to/S001_R2.fq.gz
S002,BATCH_1,/path/to/S002_R1.fq.gz,/path/to/S002_R2.fq.gz
```

You can auto-generate the sample sheet from a directory structured as `BATCH/SAMPLE_ID/`:

```bash
bash make_sample_sheet.sh /path/to/data/ samples.csv
```

Supported read filename patterns: `_R1.fq.gz`, `_R1.fastq.gz`, `_1.fq.gz`, `_1.fastq.gz` (and R2 equivalents).

### Running GTDB-Tk on existing MAGs

If you already have MAGs, use the `GTDBTK_ONLY` entry point with a bins sheet CSV (`sample_id`, `batch`, `bins_dir`):

```csv
sample_id,batch,bins_dir
S001,BATCH_1,/path/to/S001/bins
S002,BATCH_1,/path/to/S002/bins
```

Auto-generate from a `BATCH/SAMPLE_ID/` directory containing `.fa` files:

```bash
bash make_bins_sheet.sh /path/to/mags/ bins.csv
```

Then run:

```bash
nextflow run main.nf -entry GTDBTK_ONLY \
  --bins_sheet bins.csv \
  --gtdbtk_db /path/to/gtdbtk_data \
  --output_dir results \
  -profile standard
```

### Database requirements

**GTDB-Tk reference data** (~85 GB for r220):
```bash
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz
# Pass at runtime: --gtdbtk_db /path/to/release220/
```
> Requires ~200 GB RAM. Use `--skip_gtdbtk` to skip.

**CheckM2 database** (~3 GB):
```bash
checkm2 database --download --path /path/to/checkm2_db
# Pass at runtime: --checkm2_db /path/to/checkm2_db/uniref100.KO.1.dmnd
# Or set: export CHECKM2DB=/path/to/checkm2_db/uniref100.KO.1.dmnd
```

**Host genome Bowtie2 index** (optional):
```bash
bowtie2-build host_genome.fa /path/to/host_index/host
# Pass the directory at runtime: --host_index /path/to/host_index/
```

### Run

```bash
# Local execution (conda), with host removal
nextflow run main.nf \
  --sample_sheet samples.csv \
  --output_dir results \
  --host_index /path/to/host_index \
  --gtdbtk_db /path/to/gtdbtk_data \
  --checkm2_db /path/to/checkm2_db/uniref100.KO.1.dmnd \
  -profile standard

# HPC with PBS scheduler
nextflow run main.nf \
  --sample_sheet samples.csv \
  --output_dir results \
  --gtdbtk_db /path/to/gtdbtk_data \
  -profile pbs \
  -resume

# Skip host removal and GTDB-Tk (binning + quality only)
nextflow run main.nf \
  --sample_sheet samples.csv \
  --output_dir results \
  --skip_gtdbtk \
  -profile pbs

# Tune binning parameters
nextflow run main.nf \
  --sample_sheet samples.csv \
  --output_dir results \
  --gtdbtk_db /path/to/gtdbtk_data \
  --metabat2_min_contig 2000 \
  --maxbin2_min_contig 2000 \
  --concoct_min_contig 1000 \
  --dastool_score_threshold 0.6 \
  -profile pbs
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sample_sheet` | required | CSV with sample_id, batch, read1, read2 |
| `--bins_sheet` | null | CSV with sample_id, batch, bins_dir (GTDBTK_ONLY entry only) |
| `--output_dir` | `./results` | Output directory |
| `--fastp_min_qual` | `20` | fastp minimum base quality |
| `--fastp_min_len` | `50` | fastp minimum read length after trimming |
| `--fastp_cpus` | `8` | CPUs for fastp |
| `--host_index` | `null` | Directory with Bowtie2 host genome index (omit to skip) |
| `--bowtie2_cpus` | `24` | CPUs for Bowtie2 |
| `--min_contig_len` | `500` | Minimum contig length kept by MEGAHIT |
| `--megahit_cpus` | `24` | CPUs for MEGAHIT |
| `--metabat2_min_contig` | `1500` | MetaBat2 minimum contig length |
| `--metabat2_cpus` | `16` | CPUs for MetaBat2 |
| `--maxbin2_min_contig` | `1500` | MaxBin2 minimum contig length |
| `--maxbin2_cpus` | `16` | CPUs for MaxBin2 |
| `--concoct_chunk_size` | `10000` | CONCOCT contig chunk size (bp) |
| `--concoct_min_contig` | `1000` | CONCOCT minimum contig length |
| `--concoct_cpus` | `16` | CPUs for CONCOCT |
| `--dastool_score_threshold` | `0.5` | DAS Tool bin score threshold |
| `--dastool_cpus` | `16` | CPUs for DAS Tool |
| `--checkm2_db` | `null` | Path to CheckM2 diamond database |
| `--checkm2_cpus` | `16` | CPUs for CheckM2 |
| `--gtdbtk_db` | required | Path to GTDB-Tk reference data (unless `--skip_gtdbtk`) |
| `--gtdbtk_skani_sketch_dir` | `null` | Path to skani sketch dir for ANI screening (skips ANI screen if not provided) |
| `--gtdbtk_cpus` | `32` | CPUs for GTDB-Tk |
| `--skip_gtdbtk` | `false` | Skip GTDB-Tk taxonomic classification |
| `--skip_summary` | `false` | Skip summary generation |
| `--queue_size` | `2` | Max concurrent jobs |

## Output Structure

```
results/
└── BATCH/
    └── SAMPLE_ID/
        ├── fastp/
        │   ├── *_R1_trimmed.fq.gz
        │   ├── *_R2_trimmed.fq.gz
        │   ├── *_fastp.json
        │   └── *_fastp.html
        ├── host_removal/              # Only if --host_index is set
        │   ├── *_clean_R1.fq.gz
        │   ├── *_clean_R2.fq.gz
        │   └── *_host_removal_summary.txt
        ├── megahit/
        │   ├── *_contigs.fa
        │   └── *_megahit.log
        ├── mapping/
        │   ├── *_sorted.bam
        │   ├── *_sorted.bam.bai
        │   └── *_mapping_summary.txt
        ├── metabat2/
        │   ├── *_metabat2_bins/       # MetaBat2 bin FASTA files
        │   ├── *_depth.txt
        │   └── *_metabat2.log
        ├── maxbin2/
        │   ├── *_maxbin2_bins/        # MaxBin2 bin FASTA files
        │   └── *_maxbin2.log
        ├── concoct/
        │   ├── *_concoct_bins/        # CONCOCT bin FASTA files
        │   └── *_concoct_output/      # CONCOCT intermediate files
        ├── dastool/
        │   ├── *_dastool_bins/        # Refined MAG bins (FASTA)
        │   ├── *_dastool_DASTool_summary.tsv
        │   └── *_dastool_DASTool_scores.tsv
        ├── checkm2/
        │   └── quality_report.tsv
        └── gtdbtk/
            ├── *_bac120.summary.tsv
            └── *_ar53.summary.tsv
results/
└── summary/
    ├── all_bins_summary.tsv
    ├── all_checkm2_quality.tsv
    ├── all_gtdbtk_bac120_summary.tsv
    └── all_gtdbtk_ar53_summary.tsv
```

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution with conda/mamba |
| `pbs` | HPC with PBS/Torque scheduler and conda |
| `modules` | HPC with module system |
| `docker` | Docker containers |
| `singularity` | Singularity containers |
| `test` | Reduced resources for testing |

## File Structure

```
MAG-pipeline/
├── main.nf                     # Main workflow
├── nextflow.config             # Configuration and profiles
└── modules/
    ├── fastp.nf                # QC & trimming
    ├── bowtie2_host.nf         # Host read removal
    ├── megahit.nf              # Assembly
    ├── bowtie2_mapping.nf      # Read mapping back to assembly
    ├── metabat2.nf             # Depth calculation and MetaBat2 binning
    ├── maxbin2.nf              # MaxBin2 binning
    ├── concoct.nf              # CONCOCT binning
    ├── dastool.nf              # Bin refinement with DAS Tool
    ├── checkm2.nf              # Bin quality assessment
    ├── gtdbtk.nf               # Taxonomic classification of MAGs
    └── summarize.nf            # Per-sample and merged summaries
```

## Citations

If you use this pipeline, please cite the following tools:

**QC & Trimming**
- **fastp**: Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*. 2018;34(17):i884–i890. https://doi.org/10.1093/bioinformatics/bty560

**Host Removal & Read Mapping**
- **Bowtie2**: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nat Methods*. 2012;9(4):357–359. https://doi.org/10.1038/nmeth.1923
- **Samtools**: Li H, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 2009;25(16):2078–2079. https://doi.org/10.1093/bioinformatics/btp352

**Assembly**
- **MEGAHIT**: Li D, et al. MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*. 2015;31(10):1674–1676. https://doi.org/10.1093/bioinformatics/btv033

**Binning**
- **MetaBat2**: Kang DD, et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*. 2019;7:e7359. https://doi.org/10.7717/peerj.7359
- **MaxBin2**: Wu YW, et al. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics*. 2016;32(4):605–607. https://doi.org/10.1093/bioinformatics/btv638
- **CONCOCT**: Alneberg J, et al. Binning metagenomic contigs by coverage and composition. *Nat Methods*. 2014;11(11):1144–1146. https://doi.org/10.1038/nmeth.3103

**Bin Refinement**
- **DAS Tool**: Sieber CMK, et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol*. 2018;3(7):836–843. https://doi.org/10.1038/s41564-018-0171-1

**Quality Assessment**
- **CheckM2**: Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods*. 2023;20(8):1203–1212. https://doi.org/10.1038/s41592-023-01940-w

**Taxonomic Classification**
- **GTDB-Tk v2**: Chaumeil PA, et al. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics*. 2022;38(23):5315–5316. https://doi.org/10.1093/bioinformatics/btac672
- **GTDB**: Parks DH, et al. A complete domain-to-species taxonomy for Bacteria and Archaea. *Nat Biotechnol*. 2022;40(7):1021–1028. https://doi.org/10.1038/s41587-020-0501-8
