# Cecilia-nf

**usearCh basEd ampliCon pIpeLine for Illumina dAta** — Nextflow DSL2 conversion of Cecilia v2.0.

> Gian M.N. Benucci, Ph.D. | benucci@msu.edu | Michigan State University

---

## Overview

Cecilia-nf processes paired-end (or single-end) Illumina amplicon reads through quality control, assembly, clustering, and taxonomy. All steps run as SLURM jobs via the Nextflow executor with Singularity containers, providing per-sample parallelism and full reproducibility.

The pipeline replaces the 18 sequential SLURM bash scripts in `code/` with a single DAG-driven workflow. The original `config.yaml` is replaced by `nextflow.config`.

---

## Pipeline DAG

```
rawdata/*.fastq.gz
        │
        ▼
 01  DECOMPRESS          per sample — handles .gz / .bz2 / plain; renames 1_1/1_2 → R1/R2
        │
        ├──────────────▶ 02  FASTQC              raw read QC (parallel branch, non-blocking)
        │
        ▼
 03  REMOVE_PHIX         USEARCH -filter_phix, per sample
        │
        ▼
 04  ASSEMBLE_READS *    USEARCH -fastq_mergepairs, per sample
        │
        ▼
 05  RENAME_READS *      collect all → sequential sample001/002/… headers; re-parallelise
        │
        ▼
 06  STRIP_PRIMERS *     Cutadapt (-g fwd -a rev_rc -n 2), per sample
        │
        ▼
 07  POOL_READS          cat all stripped FASTQs → pooled.fastq; seqtk subset_500.fasta
        │
        ▼
 08  EE_STATS            USEARCH -fastq_eestats2; optional stripleft trim
        │
        ├──────────────▶ 08b TEST_LENGTH *        150–300 bp sweep → testLength.R plots
        │
        ▼  (× N lengths)
 09  FILTER_MAXEE        USEARCH -fastq_filter per truncation length + FastQC; optional output
        │
        ▼
 10  DEREPLICATE         USEARCH -fastx_uniques; SWARM 0-diff linearisation (if cluster_swarm)
        │
        ├──────────────▶ 11  CLUSTER_UPARSE *    USEARCH -cluster_otus + OTU table
        ├──────────────▶ 12  CLUSTER_UNOISE *    USEARCH -unoise3 → ASVs + ASV table
        │                        │
        │                        └──────────────▶ 13  CLUSTER_ASV_TO_OTU *  -cluster_smallmem
        ├──────────────▶ 14  CLOSED_REF_OTU *    USEARCH -closed_ref vs closedRef_db
        ├──────────────▶ 15  CLUSTER_SWARM *     SWARM d=1 fastidious → seeds → OTU table
        │
        ▼  (all cluster FASTAs)
 17  BUILD_SINTAX_DB *   USEARCH -makeudb_usearch (once, shared)
        │
        ▼  (× cluster FASTA)
 17  TAXONOMY_SINTAX *   USEARCH -sintax per FASTA, parallel
        │
        ▼
 18  COLLECT_RESULTS     tar.gz archive of all outputs, dated
```

Steps marked `*` are optional and controlled by flags in `nextflow.config`.

---

## Repository structure

```
cecilia-nf/
├── main.nf                  # Entry point — prints run summary, calls CECILIA workflow
├── nextflow.config          # All user parameters, SLURM resources, Singularity config
├── run_cecilia.sb           # SLURM launcher for the Nextflow head process
├── workflows/
│   ├── cecilia.nf           # Main DAG wiring all modules together
│   └── test_length.nf       # Optional length-sweep subworkflow
├── modules/                 # One .nf file per process (18 total)
│   ├── decompress.nf
│   ├── fastqc.nf
│   ├── remove_phix.nf
│   ├── assemble_reads.nf
│   ├── rename_reads.nf
│   ├── strip_primers.nf
│   ├── pool_reads.nf
│   ├── ee_stats.nf
│   ├── test_length.nf
│   ├── filter_maxee.nf
│   ├── dereplicate.nf
│   ├── cluster_uparse.nf
│   ├── cluster_unoise.nf
│   ├── cluster_asv_to_otu.nf
│   ├── closed_ref_otu.nf
│   ├── cluster_swarm.nf
│   ├── build_sintax_db.nf
│   ├── taxonomy_sintax.nf
│   └── collect_results.nf
└── bin/                     # Helper scripts auto-added to PATH in every container
    ├── renameSeq.sh
    ├── compareSeqQual.sh
    └── testLength.R
```

---

## Prerequisites

| Requirement | Notes |
|---|---|
| Nextflow ≥ 24 | `conda activate nextflow` on the cluster |
| Singularity / Apptainer | Available on SLURM compute nodes |
| USEARCH 11 | Proprietary binary, bind-mounted from host (not containerised) |
| SLURM | Default executor; `local` and `test` profiles also available |

All other tools (FastQC, Cutadapt, seqtk, SWARM, R/tidyverse) are pulled automatically as Singularity images on first run and cached in `singularity.cacheDir`.

---

## Quick start

```bash
# 1. Place paired FASTQ files (*.fastq.gz) in rawdata/
#    Expected naming: SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz

# 2. Edit nextflow.config — update paths that are site-specific:
#    params.usearch    — path to USEARCH binary
#    params.sintax_db  — path to SILVA / UNITE FASTA for taxonomy
#    params.closedRef_db — path to closed-reference FASTA (if used)

# 3. Submit to SLURM from a development node
sbatch cecilia-nf/run_cecilia.sb

# 4. Resume a previous run after a failure
sbatch cecilia-nf/run_cecilia.sb -resume <session-id>
```

The session ID is shown in the Nextflow run log and in `.nextflow/history`.

---

## Configuration

All user-facing parameters are in the `params {}` block of `nextflow.config`. Any parameter can be overridden at the command line without editing the file:

```bash
sbatch run_cecilia.sb \
  --dna_marker ITS \
  --fwd_primer CTTGGTCATTTAGAGGAAGTAA \
  --rev_primer TCCTCCGCTTATTGATATGC \
  --rev_primer_rc GCATATCAATAAGCGGAGGA \
  --max_eerr 0.5 \
  --cluster_swarm false \
  --sintax_db /path/to/UNITE.fasta
```

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `rawdata` | `../rawdata` | Directory containing input FASTQ files |
| `outdir` | `../outputs` | Results destination |
| `usearch` | *(site path)* | USEARCH binary (host path, bind-mounted) |
| `sintax_db` | *(site path)* | FASTA database for SINTAX taxonomy |
| `closedRef_db` | *(site path)* | FASTA database for closed-reference OTUs |
| `dna_marker` | `16S` | Label used in output file names |
| `paired` | `true` | `true` = R1+R2 paired-end; `false` = single-end |
| `assemble` | `true` | Merge paired reads with USEARCH |
| `rename` | `true` | Rename reads to sequential sample001/002/… |
| `primers` | `true` | Strip primers with Cutadapt |
| `stripleft` | `false` | Trim fixed N bp from read left (set `stripleft_bp`) |
| `test_length` | `false` | Run 150–300 bp length-sweep subworkflow |
| `filter_lengths` | `[200, 225, 250]` | Truncation lengths for FILTER_MAXEE |
| `user_length` | `290` | Extra truncation length |
| `max_eerr` | `1.0` | Maximum expected errors for USEARCH filter |
| `cluster_otu` | `true` | UPARSE OTU clustering |
| `cluster_asv` | `true` | UNOISE3 ASV denoising |
| `cluster_asv_to_otu` | `true` | Cluster ASVs into OTUs at `clustering_id` |
| `cluster_swarm` | `true` | SWARM d=1 fastidious clustering |
| `closed_ref_otu` | `true` | Closed-reference OTU picking |
| `sintax_taxonomy` | `true` | SINTAX taxonomy assignment |
| `otu_min_size` | `1` | Minimum OTU size (UPARSE) |
| `zotu_min_size` | `8` | Minimum ZOTU size (UNOISE3) |
| `clustering_id` | `0.97` | Identity threshold for ASV→OTU clustering |
| `sintax_cutoff` | `0.8` | SINTAX confidence cutoff |
| `fwd_primer` | `GTGCCAGCMGCCGCGGTAA` | 515F (Parada/Caporaso) |
| `rev_primer` | `GGACTACHVGGGTWTCTAAT` | 806R (Apprill/Caporaso) |
| `slurm_account` | `glbrc` | SLURM billing account |

---

## Outputs

```
outputs/
├── 01_decompressed/          raw decompressed FASTQs
├── 02_rawQuality_fqc/        FastQC HTML + ZIP reports
├── 03_removePhix_usearch/    PhiX-filtered FASTQs
├── 04_readAssembly_usearch/  assembled (merged) FASTQs
├── 05_readRenamed_bash/      renamed FASTQs + sample mapping file
├── 06_stripPrimers_cutadapt/ primer-stripped FASTQs
├── 07_poolReads/             pooled.fastq + subset_500.fasta
├── 08_eeStats_usearch/       eestats2 tables + compareSeqQual output
├── 09_filterMaxee_usearch/   filtered FASTQs + FastQC per length
├── 10_dereplicate_usearch/   uniques.fasta per length
├── 11_clusterOTU_usearch/    OTU FASTA + OTU table (UPARSE)
├── 12_clusterASV_usearch/    ASV FASTA + ASV table (UNOISE3)
├── 13_taxonomy_sintax/       SINTAX taxonomy per cluster set
├── 14_closedRef_usearch/     closed-reference OTU table
├── 15_clusterSwarm/          SWARM FASTA + OTU table
├── stats/                    per-sample read counts at each step
├── cecilia_report.html       Nextflow execution report
├── cecilia_timeline.html     Nextflow timeline
└── cecilia_dag.html          Nextflow workflow DAG
```

---

## Containers

| Process(es) | Image |
|---|---|
| DECOMPRESS, REMOVE_PHIX, ASSEMBLE_READS, RENAME_READS, EE_STATS, all clustering, BUILD_SINTAX_DB, TAXONOMY_SINTAX, COLLECT_RESULTS | `docker://ubuntu:22.04` + USEARCH bind-mount |
| FASTQC, FILTER_MAXEE | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| STRIP_PRIMERS | `quay.io/biocontainers/cutadapt:5.2--py312h247cb63_1` |
| POOL_READS | `quay.io/biocontainers/seqtk:1.4--he4a0461_1` |
| DEREPLICATE, CLUSTER_SWARM | `quay.io/biocontainers/swarm:3.1.5--h9f5acd7_0` |
| TEST_LENGTH | `docker://rocker/tidyverse:4.3.1` |

Images are pulled once and cached in `singularity.cacheDir` (see `nextflow.config`).

---

## Implementation notes

- **USEARCH is not containerised.** It is a proprietary binary bind-mounted from the host filesystem into every Singularity container via `singularity.runOptions`. The path is set in `params.usearch`.
- **RENAME_READS is a collect step.** All assembled FASTQs are gathered before renaming so that sequential `sample001/002/…` numbering is assigned correctly. Outputs are re-parallelised with `.flatten().map()`.
- **FILTER_MAXEE uses `optional: true`.** If no reads survive filtering at a given length, that length produces no output and is silently skipped — downstream steps do not fail.
- **SINTAX uses `.combine()`** to fan out taxonomy jobs: one SLURM job per cluster FASTA × one shared `.udb` file built once by BUILD_SINTAX_DB.
- **FastQC on parallel filesystems.** JVM bus errors (SIGBUS) on GPFS/BeeGFS are fixed by setting `_JAVA_OPTIONS="-Djava.io.tmpdir=/tmp -Xss512k"` and using `--dir /tmp` in the FastQC call. The FASTQC process also has `errorStrategy = 'retry'; maxRetries = 1`.
- **CONSTAX2 is intentionally excluded.** Use SINTAX taxonomy output if needed.
