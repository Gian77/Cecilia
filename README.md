<p align="center">
  <img alt="Cecilia-nf logo" src=".logo/cecilia_logo.png" width="100%">
</p>

A Nextflow DSL2 reimplementation of Cecilia v2.0.

> Gian M.N. Benucci, Ph.D. | benucci@msu.edu | Michigan State University

---

## Overview

CECILIA processes paired-end (or single-end) Illumina amplicon reads through quality control, assembly, clustering, and taxonomy. All steps run as SLURM jobs via the Nextflow executor with Singularity containers, providing per-sample parallelism and full reproducibility.

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
CECILIAv2/
├── cecilia.nf               # Entry point and main DAG wiring all modules together
├── nextflow.config          # All user parameters, SLURM resources, Singularity config
├── run_cecilia.sb           # SLURM launcher for the Nextflow head process
├── lib/
│   └── Utils.groovy         # Shared helpers (param coercion, USEARCH fetch)
├── modules/                 # One .nf file per process (20 total)
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
│   ├── summarize_counts.nf
│   └── collect_results.nf
└── bin/                     # Helper scripts auto-added to PATH in every container
    ├── renameSeq.sh
    ├── compareSeqQual.sh
    ├── summarizeCounts.sh
    └── testLength.R
```

---

## Prerequisites

| Requirement | Notes |
|---|---|
| Nextflow ≥ 24 | `conda activate nextflow` on the cluster |
| Singularity / Apptainer | Available on SLURM compute nodes |
| USEARCH 11.0.667 | Public domain (CC0-1.0); downloaded automatically on first run, bind-mounted from host (not containerised). Version is pinned deliberately — v12 was evaluated and rejected, see `nextflow.config` |
| SLURM | Default executor; `local` and `test` profiles also available |

All other tools (FastQC, Cutadapt, seqtk, SWARM, R/tidyverse) are pulled automatically as Singularity images on first run and cached in `singularity.cacheDir`.

---

## Quick start

```bash
# 1. Place paired FASTQ files (*.fastq.gz) in rawdata/
#    Expected naming: SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz

# 2. Edit nextflow.config — update paths that are site-specific:
#    params.sintax_db  — path to SILVA / UNITE FASTA for taxonomy
#    params.closedRef_db — path to closed-reference FASTA (if used)
#    USEARCH needs no setup: it is downloaded automatically on first run.

# 3a. Run interactively from a development node (inside tmux)
cd CECILIAv2
module purge && module load Miniforge3 && conda activate nextflow
nextflow run cecilia.nf

# 3b. ...or submit the head process itself to SLURM
sbatch CECILIAv2/run_cecilia.sb

# 4. Resume a previous run after a failure (either form)
nextflow run cecilia.nf -resume <session-id>
sbatch CECILIAv2/run_cecilia.sb -resume <session-id>
```

The session ID is shown in the Nextflow run log and in `.nextflow/history`. See
[Configuration](#configuration) for the difference between the two launch modes and
for passing parameters on the command line.

---

## Configuration

All user-facing parameters are in the `params {}` block of `nextflow.config`. Any parameter can be overridden at the command line without editing the file. The same flags work whether you launch from the terminal or through SLURM.

### Running from the terminal

Run the Nextflow head process interactively on a **development node** — never a login/gateway node. Individual pipeline steps are still submitted to SLURM (`process.executor = 'slurm'` is the default), so the terminal process only orchestrates and stays lightweight:

```bash
ssh dev-amd20                     # any ICER development node
cd /path/to/CECILIAv2

module purge
module load Miniforge3
conda activate nextflow

nextflow run cecilia.nf \
  --dna_marker ITS \
  --fwd_primer CTTGGTCATTTAGAGGAAGTAA \
  --rev_primer TCCTCCGCTTATTGATATGC \
  --rev_primer_rc GCATATCAATAAGCGGAGGA \
  --max_eerr 0.5 \
  --cluster_swarm false \
  --sintax_db /path/to/UNITE.fasta
```

The head process has to stay alive for the entire run, so start it inside `tmux` (or `screen`) — otherwise it is killed the moment your SSH session drops:

```bash
tmux new -s cecilia
# ... launch nextflow, then detach with Ctrl-b then d
tmux attach -t cecilia            # reattach later, from any login node
```

Print the full parameter list and exit, without launching anything:

```bash
nextflow run cecilia.nf --help
```

Useful Nextflow flags (single dash — distinct from pipeline `--params`):

| Flag | Effect |
|---|---|
| `-resume` | reuse cached results from the previous run |
| `-profile local` | run every process on the current node instead of submitting to SLURM |
| `-profile test` | as `local`, with `test_length` and `cluster_swarm` forced off |
| `-ansi-log false` | plain scrolling log, better when piping to `tee` or a file |
| `-with-report` | write an execution report (already enabled by default in `nextflow.config`) |

Note that `local` and `test` also disable Singularity, so every tool (FastQC, Cutadapt, seqtk, SWARM, R) must already be on your `PATH`. They are meant for quick wiring and syntax checks, not for real runs.

### Running as a batch job

For long runs, submit the head process itself to SLURM so it survives disconnection with no `tmux` needed. Parameters are passed straight through to `nextflow run`:

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

Progress goes to `cecilia_nf-<jobid>.log` in the launch directory. The walltime in `run_cecilia.sb` (`--time=03:59:00`) caps the head process, not the individual steps — raise it for large runs.

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `rawdata` | `../rawdata` | Directory containing input FASTQ files |
| `outdir` | `../outputs` | Results destination |
| `usearch` | `.tools/usearch11.0.667_i86linux64` | USEARCH binary; auto-downloaded on first run |
| `usearch_fetch` | `true` | Allow the download; `false` = fail if the binary is absent |
| `sintax_db` | *(site path)* | FASTA database for SINTAX taxonomy |
| `closedRef_db` | *(site path)* | FASTA database for closed-reference OTUs |
| `dna_marker` | `16S` | Label used in output file names |
| `paired` | `true` | `true` = R1+R2 paired-end; `false` = single-end |
| `assemble` | `true` | Merge paired reads with USEARCH |
| `rename` | `true` | Rename reads to sequential sample001/002/… |
| `primers` | `true` | Strip primers with Cutadapt |
| `stripleft` | `false` | Trim fixed N bp from read left (set `stripleft_bp`) |
| `test_length` | `false` | Run the 150–300 bp length sweep |
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

### Primer pairs

Common primer pairs, carried over from the legacy `Cecilia/config.yaml`, so you don't have to go hunting for sequences. Set `fwd_primer` / `rev_primer` / `rev_primer_rc` on the command line to switch marker/pair.

| Marker | Pair | `fwd_primer` | `rev_primer` | `rev_primer_rc` | Reference |
|---|---|---|---|---|---|
| 16S | 515F / 806R | `GTGCCAGCMGCCGCGGTAA` | `GGACTACHVGGGTWTCTAAT` | `ATTAGAWACCCBDGTAGTCC` | Parada et al. 2016, modified by Caporaso; Apprill et al. 2015, modified by Caporaso ([EMP protocol](https://earthmicrobiome.org/protocols-and-standards/16s/)) — **pipeline default** |
| 16S | 27F / 1492R | `AGAGTTTGATCMTGGCTCAG` | `GGTTACCTTGTTACGACTT` | `AAGTCGTAACAAGGTAACC` | Heuer et al., *Appl Environ Microbiol* 1997;63:3233–41 |
| 16S | 341F / 806R | `CCTACGGGNGGCWGCAG` (or `CCTACGGGAGGCAGCAG`) | `GGACTACHVGGGTWTCTAAT` | `ATTAGAWACCCBDGTAGTCC` | Albertsen et al., *PLoS ONE* 2015;10:e0132783 |
| ITS | ITS1F / ITS4 | `CTTGGTCATTTAGAGGAAGTAA` | `TCCTCCGCTTATTGATATGC` | `GCATATCAATAAGCGGAGGA` | Gardes & Bruns 1993; White et al. 1990 |
| ITS | ITS1F / ITS2 | `CTTGGTCATTTAGAGGAAGTAA` | `GCTGCGTTCTTCATCGATGC` | `GCATCGATGAAGAACGCAGC` | Gardes & Bruns 1993; White et al. 1990 |
| ITS | 5.8S-Fun / ITS4-Fun | `AACTTTYRRCAAYGGATCWCT` | `AGCCTCCGCTTATTGATATGCTTAART` | `AYTTAAGCATATCAATAAGCGGAGGCT` | Taylor et al., *Appl Environ Microbiol* 2016 |
| ITS | ITS5-1737F / ITS2-2043R | `GGAAGTAAAAGTCGTAACAAGG` | `GCTGCGTTCTTCATCGATGC` | `GCATCGATGAAGAACGCAGC` | Novogene ITS primer set |
| AMF (18S) | SSUmCf / LSUmBr | `TATYGYTCTTNAACGAGGAATC` | `AACACTCGCAYAYATGYTAGA` | `TCTARCATRTRTGCGAGTGTT` | Krüger et al. 2009 |

---

## References

Full citations for the primer pairs in the table above:

- Parada, A.E., Needham, D.M. & Fuhrman, J.A. (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. *Environmental Microbiology*, 18(5), 1403–1414. https://doi.org/10.1111/1462-2920.13023
- Apprill, A., McNally, S., Parsons, R. & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. *Aquatic Microbial Ecology*, 75, 129–137. https://doi.org/10.3354/ame01753
- Heuer, H., Krsek, M., Baker, P., Smalla, K. & Wellington, E.M.H. (1997). Analysis of actinomycete communities by specific amplification of genes encoding 16S rRNA and gel-electrophoretic separation in denaturing gradients. *Applied and Environmental Microbiology*, 63(8), 3233–3241.
- Albertsen, M., Karst, S.M., Ziegler, A.S., Kirkegaard, R.H. & Nielsen, P.H. (2015). Back to basics – the influence of DNA extraction and primer choice on phylogenetic analysis of activated sludge communities. *PLoS ONE*, 10(7), e0132783. https://doi.org/10.1371/journal.pone.0132783
- Gardes, M. & Bruns, T.D. (1993). ITS primers with enhanced specificity for basidiomycetes—application to the identification of mycorrhizae and rusts. *Molecular Ecology*, 2(2), 113–118. https://doi.org/10.1111/j.1365-294X.1993.tb00005.x
- White, T.J., Bruns, T., Lee, S. & Taylor, J. (1990). Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics. In Innis, M.A., Gelfand, D.H., Sninsky, J.J. & White, T.J. (Eds.), *PCR Protocols: A Guide to Methods and Applications* (pp. 315–322). Academic Press.
- Taylor, D.L., Walters, W.A., Lennon, N.J., Bochicchio, J., Krohn, A., Caporaso, J.G. & Pennanen, T. (2016). Accurate estimation of fungal diversity and abundance through improved lineage-specific primers optimized for Illumina amplicon sequencing. *Applied and Environmental Microbiology*, 82(24), 7217–7226. https://doi.org/10.1128/AEM.02576-16
- Krüger, M., Stockinger, H., Krüger, C. & Schüßler, A. (2009). DNA-based species level detection of Glomeromycota: one PCR primer set for all arbuscular mycorrhizal fungi. *New Phytologist*, 183(1), 212–223. https://doi.org/10.1111/j.1469-8137.2009.02835.x
- ITS5-1737F / ITS2-2043R: proprietary Novogene ITS primer set, no peer-reviewed citation available.

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

- **USEARCH is fetched, not containerised.** USEARCH v11.0.667 was donated to the public domain under CC0-1.0 in 2024 ([`rcedgar/usearch_old_binaries`](https://github.com/rcedgar/usearch_old_binaries)), so the pipeline downloads it into `.tools/` on first run and bind-mounts it into every Singularity container via `singularity.runOptions`. Nothing for the user to obtain by hand. The URL is pinned to a commit rather than a branch, and the download is SHA-256 verified before it is moved into place, so an interrupted fetch cannot leave a truncated binary behind. Set `--usearch /path/to/binary` to use a local copy instead — an existing file is never overwritten, and a checksum mismatch there only warns.
- **The fetch runs in the entry workflow, not in a process.** Every USEARCH module interpolates `${params.usearch}` straight into its script, and Nextflow params cannot be reassigned at runtime — `params.usearch = …` from a workflow body is silently ignored (verified on 26.04.3). The binary therefore has to exist before the first task is submitted, which is what `Utils.ensureUsearch()` guarantees. Passing it through a channel instead would mean rewriting the input block of all fifteen USEARCH modules and invalidating the entire `-resume` cache.
- **v11, deliberately, not v12.** USEARCH v12 is open source (GPL-3.0) but drops five commands this pipeline needs (`-filter_phix`, `-fastq_eestats2`, `-sortbylength`, `-fastx_getseqs`, `-fastx_info`), and its `-closed_ref` segfaults on input that v11 processes normally. Since v11 is now public domain, v12 offers no licensing advantage either.
- **RENAME_READS is a collect step.** All assembled FASTQs are gathered before renaming so that sequential `sample001/002/…` numbering is assigned correctly. Outputs are re-parallelised with `.flatten().map()`.
- **FILTER_MAXEE uses `optional: true`.** If no reads survive filtering at a given length, that length produces no output and is silently skipped — downstream steps do not fail.
- **SINTAX uses `.combine()`** to fan out taxonomy jobs: one SLURM job per cluster FASTA × one shared `.udb` file built once by BUILD_SINTAX_DB.
- **FastQC on parallel filesystems.** JVM bus errors (SIGBUS) on GPFS/BeeGFS are fixed by setting `_JAVA_OPTIONS="-Djava.io.tmpdir=/tmp -Xss512k"` and using `--dir /tmp` in the FastQC call. The FASTQC process also has `errorStrategy = 'retry'; maxRetries = 1`.
- **CONSTAX2 is intentionally excluded.** Use SINTAX taxonomy output if needed.
