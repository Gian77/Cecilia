#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// =============================================================================
// Cecilia-nf  —  cecilia.nf
// usearCh basEd ampliCon pIpeLine for Illumina dAta
// Gian M.N. Benucci, Ph.D.  |  benucci@msu.edu
// =============================================================================

include { DECOMPRESS        } from './modules/decompress'
include { FASTQC            } from './modules/fastqc'
include { REMOVE_PHIX       } from './modules/remove_phix'
include { ASSEMBLE_READS    } from './modules/assemble_reads'
include { RENAME_READS      } from './modules/rename_reads'
include { STRIP_PRIMERS     } from './modules/strip_primers'
include { POOL_READS        } from './modules/pool_reads'
include { EE_STATS          } from './modules/ee_stats'
include { FILTER_MAXEE      } from './modules/filter_maxee'
include { DEREPLICATE       } from './modules/dereplicate'
include { CLUSTER_UPARSE    } from './modules/cluster_uparse'
include { CLUSTER_UNOISE    } from './modules/cluster_unoise'
include { CLUSTER_ASV_TO_OTU } from './modules/cluster_asv_to_otu'
include { CLOSED_REF_OTU    } from './modules/closed_ref_otu'
include { CLUSTER_SWARM     } from './modules/cluster_swarm'
include { BUILD_SINTAX_DB   } from './modules/build_sintax_db'
include { TAXONOMY_SINTAX   } from './modules/taxonomy_sintax'
include { SUMMARIZE_COUNTS  } from './modules/summarize_counts'
include { COLLECT_RESULTS   } from './modules/collect_results'
include { TEST_LENGTH       } from './modules/test_length'

// =============================================================================
// MAIN CECILIA WORKFLOW
// =============================================================================
workflow CECILIA {

    // -------------------------------------------------------------------------
    // Boolean params arrive as Strings when overridden on the CLI (e.g.
    // `--paired false`) — coerce once here and use these locals everywhere
    // below instead of the raw params.* to avoid Groovy string-truthiness bugs.
    // -------------------------------------------------------------------------
    def paired             = Utils.asBool(params.paired)
    def assemble           = Utils.asBool(params.assemble)
    def rename              = Utils.asBool(params.rename)
    def stripleft            = Utils.asBool(params.stripleft)
    def test_length          = Utils.asBool(params.test_length)
    def cluster_otu           = Utils.asBool(params.cluster_otu)
    def cluster_asv           = Utils.asBool(params.cluster_asv)
    def cluster_asv_to_otu    = Utils.asBool(params.cluster_asv_to_otu)
    def closed_ref_otu        = Utils.asBool(params.closed_ref_otu)
    def cluster_swarm         = Utils.asBool(params.cluster_swarm)
    def sintax_taxonomy       = Utils.asBool(params.sintax_taxonomy)
    def summary_table         = Utils.asBool(params.summary_table)

    // -------------------------------------------------------------------------
    // INPUT CHANNEL
    // Paired-end: fromFilePairs matches *_R1_* / *_R2_* or *_1.fastq / *_2.fastq
    // Single-end: fromPath with a per-sample tuple
    // -------------------------------------------------------------------------
    if (paired) {
        ch_raw = Channel
            .fromFilePairs("${params.rawdata}/*_{R1,R2}*.fastq{,.gz,.bz2}", checkIfExists: true)
    } else {
        ch_raw = Channel
            .fromPath("${params.rawdata}/*.fastq{,.gz,.bz2}", checkIfExists: true)
            .map { f -> [ f.simpleName, [f] ] }
    }

    // -------------------------------------------------------------------------
    // STEP 01  —  Decompress  (per sample, parallel)
    // -------------------------------------------------------------------------
    DECOMPRESS(ch_raw)

    // -------------------------------------------------------------------------
    // STEP 02  —  FastQC on raw reads  (parallel branch, does not block DAG)
    // R1 and R2 are analysed independently: one SLURM job per direction,
    // each receiving all samples for that direction pooled together.
    // groupTuple() collects all per-sample files sharing the same group key.
    // -------------------------------------------------------------------------
    if (paired) {
        ch_fqc_r1 = DECOMPRESS.out.reads.map { id, reads -> tuple('R1', Utils.asList(reads)[0]) }.groupTuple()
        ch_fqc_r2 = DECOMPRESS.out.reads.map { id, reads -> tuple('R2', Utils.asList(reads)[1]) }.groupTuple()
        FASTQC(ch_fqc_r1.mix(ch_fqc_r2))
    } else {
        FASTQC(DECOMPRESS.out.reads.map { id, reads -> tuple('all', Utils.asList(reads)[0]) }.groupTuple())
    }

    // -------------------------------------------------------------------------
    // STEP 03  —  PhiX removal  (per sample)
    // Counts are collected across all samples into one summary file in stats/.
    // -------------------------------------------------------------------------
    REMOVE_PHIX(DECOMPRESS.out.reads)
    // Keep a handle on the merged counts file: SUMMARIZE_COUNTS consumes it.
    // Assignment only — the collectFile behaviour is unchanged.
    ch_nophix_counts = REMOVE_PHIX.out.counts
        .collectFile(name: 'all.nophix.counts', storeDir: "${params.outdir}/stats", newLine: false)

    // -------------------------------------------------------------------------
    // STEP 04  —  Paired-read assembly  (optional, per sample)
    // -------------------------------------------------------------------------
    if (assemble) {
        ASSEMBLE_READS(REMOVE_PHIX.out.reads)
        ch_assembled_counts = ASSEMBLE_READS.out.counts
            .collectFile(name: 'all.assembled.counts', storeDir: "${params.outdir}/stats", newLine: false)
        ch_after_assemble = ASSEMBLE_READS.out.reads
    } else {
        ch_assembled_counts = Channel.empty()
        ch_after_assemble = REMOVE_PHIX.out.reads
    }

    // -------------------------------------------------------------------------
    // STEP 05  —  Read renaming  (optional)
    // RENAME_READS is a collect step: all samples go in together so
    // sequential numbering (sample001, sample002, ...) is assigned correctly.
    // Outputs are re-parallelised with flatten().map() for downstream steps.
    // -------------------------------------------------------------------------
    if (rename) {
        RENAME_READS(ch_after_assemble.map { id, f -> f }.collect())
        ch_after_rename = RENAME_READS.out.reads
            .flatten()
            .map { f -> [ f.baseName, f ] }
    } else {
        ch_after_rename = ch_after_assemble
    }

    // -------------------------------------------------------------------------
    // STEP 06  —  Primer/adapter stripping  (per sample)
    // then pool all stripped FASTQs into one pooled.fastq
    // -------------------------------------------------------------------------
    STRIP_PRIMERS(ch_after_rename)

    POOL_READS(STRIP_PRIMERS.out.fastq.map { id, f -> f }.collect())

    // Fork the pooled channel: it is reused much later by OTU table generation
    ch_pooled = POOL_READS.out.pooled_fastq

    // -------------------------------------------------------------------------
    // STEP 07  —  EE stats  (pooled)
    // Optionally strip left bases; output may be trimmed.fastq or pooled.fastq
    // -------------------------------------------------------------------------
    EE_STATS(ch_pooled)

    // Downstream source for filtering: trimmed.fastq if stripleft, else pooled.fastq
    // EE_STATS.out.trimmed is optional; using ch_pooled when stripleft=false avoids
    // an empty-channel problem downstream.
    ch_source = stripleft ? EE_STATS.out.trimmed : ch_pooled

    // -------------------------------------------------------------------------
    // STEP 08  —  Test read length  (optional)
    // -------------------------------------------------------------------------
    if (test_length) {
        TEST_LENGTH(ch_source)
    }

    // -------------------------------------------------------------------------
    // STEP 09  —  Filter by max expected errors  (per length, parallel)
    // All lengths (fixed list + user_length) deduplicated and run concurrently.
    // FILTER_MAXEE declares its tuple output optional: processes that produce
    // no reads simply emit nothing — no downstream failure.
    // -------------------------------------------------------------------------
    def all_lengths = (params.filter_lengths + [params.user_length]).unique()
    ch_lengths = Channel.fromList(all_lengths)

    FILTER_MAXEE(ch_lengths.combine(ch_source))

    // -------------------------------------------------------------------------
    // STEP 10  —  Dereplicate  (per filtered file, parallel)
    // Emits both uniques.fasta (for UPARSE/UNOISE/closed-ref) and
    // *_linear.fasta (for SWARM)
    // -------------------------------------------------------------------------
    // Only non-empty filter outputs reach DEREPLICATE (optional: true above)
    DEREPLICATE(FILTER_MAXEE.out.filtered)

    // Attach the single pooled.fastq to each (len, uniques.fasta) pair
    // so OTU-table steps have everything they need
    ch_uniques_pool = DEREPLICATE.out.uniques.combine(ch_pooled)

    // -------------------------------------------------------------------------
    // STEPS 11–15  —  Clustering  (all optional, all parallel)
    // -------------------------------------------------------------------------
    ch_cluster_fastas = Channel.empty()

    if (cluster_otu) {
        CLUSTER_UPARSE(ch_uniques_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_UPARSE.out.fasta)
    }

    if (cluster_asv) {
        CLUSTER_UNOISE(ch_uniques_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_UNOISE.out.asv)

        if (cluster_asv_to_otu) {
            CLUSTER_ASV_TO_OTU(
                CLUSTER_UNOISE.out.asv
                    .map { len, fasta -> tuple(len, fasta) }
                    .combine(ch_pooled)
            )
            ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_ASV_TO_OTU.out.fasta)
        }
    }

    if (closed_ref_otu) {
        // Needs pooled.fastq as well: per-sample counts come from the labelled
        // pooled reads, not from the cross-sample dereplicated uniques.
        CLOSED_REF_OTU(ch_uniques_pool)
    }

    if (cluster_swarm) {
        ch_linear_pool = DEREPLICATE.out.linear.combine(ch_pooled)
        CLUSTER_SWARM(ch_linear_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_SWARM.out.fasta)
    }

    // -------------------------------------------------------------------------
    // STEP 17  —  SINTAX taxonomy  (optional)
    // BUILD_SINTAX_DB runs once; TAXONOMY_SINTAX fans out over all cluster
    // FASTAs via .combine() — one SLURM job per FASTA, all parallel.
    // -------------------------------------------------------------------------
    if (sintax_taxonomy) {
        BUILD_SINTAX_DB(file(params.sintax_db))
        TAXONOMY_SINTAX(
            ch_cluster_fastas.combine(BUILD_SINTAX_DB.out.udb)
        )
    }

    // -------------------------------------------------------------------------
    // STEP 18  —  Collect and package results
    // Assemble two flat channels:
    //   ch_misc    — QC, EE stats, subset FASTA, file mapping, testLength
    //   ch_cluster — all cluster FASTAs, OTU tables, taxonomy outputs
    // Both are .collect()ed so COLLECT_RESULTS waits for the entire pipeline.
    // -------------------------------------------------------------------------

    // Misc results (always present)
    ch_misc = Channel.empty()
        .mix(FASTQC.out.html)
        .mix(POOL_READS.out.subset_fasta)
        .mix(EE_STATS.out.eestats)

    if (rename) {
        ch_misc = ch_misc.mix(RENAME_READS.out.mapping)
    }
    if (test_length) {
        ch_misc = ch_misc
            .mix(TEST_LENGTH.out.results)
            .mix(TEST_LENGTH.out.plots)
    }
    ch_misc = ch_misc.mix(FILTER_MAXEE.out.fastqc_html)

    // Cluster results (map tuples to plain paths)
    ch_cluster = Channel.empty()
        .mix(DEREPLICATE.out.uniques.map { len, f -> f })

    if (cluster_otu) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_UPARSE.out.fasta.map { len, f -> f })
            .mix(CLUSTER_UPARSE.out.table)
    }
    if (cluster_asv) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_UNOISE.out.asv.map { len, f -> f })
            .mix(CLUSTER_UNOISE.out.table)
    }
    if (cluster_asv && cluster_asv_to_otu) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_ASV_TO_OTU.out.fasta.map { len, f -> f })
            .mix(CLUSTER_ASV_TO_OTU.out.table)
    }
    if (closed_ref_otu) {
        ch_cluster = ch_cluster
            .mix(CLOSED_REF_OTU.out.table.map { len, f -> f })
    }
    if (cluster_swarm) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_SWARM.out.fasta.map { len, f -> f })
            .mix(CLUSTER_SWARM.out.table)
    }
    if (sintax_taxonomy) {
        ch_cluster = ch_cluster
            .mix(TAXONOMY_SINTAX.out.sintax.map { len, f -> f })
    }

    // -------------------------------------------------------------------------
    // STEP 19  —  Per-sample summary tables  (optional)
    // Joins counts that upstream steps already emit into read_tracking.tsv and
    // otu_summary.tsv. Everything optional is funnelled through one collected
    // channel, so a stage that did not run simply stages no file.
    // -------------------------------------------------------------------------
    if (summary_table) {
        ch_tables = Channel.empty()
        if (cluster_otu) {
            ch_tables = ch_tables.mix(CLUSTER_UPARSE.out.table)
        }
        if (cluster_asv) {
            ch_tables = ch_tables.mix(CLUSTER_UNOISE.out.table)
        }
        if (cluster_asv && cluster_asv_to_otu) {
            ch_tables = ch_tables.mix(CLUSTER_ASV_TO_OTU.out.table)
        }
        if (closed_ref_otu) {
            ch_tables = ch_tables.mix(CLOSED_REF_OTU.out.table.map { len, f -> f })
        }
        if (cluster_swarm) {
            ch_tables = ch_tables.mix(CLUSTER_SWARM.out.table)
        }

        ch_summary_aux = Channel.empty()
            .mix(ch_nophix_counts)
            .mix(ch_assembled_counts)
            .mix(POOL_READS.out.counts)
            .mix(FILTER_MAXEE.out.filtered.map { len, f -> f })
            .mix(ch_tables)

        if (rename) {
            ch_summary_aux = ch_summary_aux.mix(RENAME_READS.out.mapping)
        }

        SUMMARIZE_COUNTS(
            file(params.rawdata),
            ch_summary_aux.collect()
        )

        // Ship the tables in the handoff archive alongside everything else
        ch_misc = ch_misc
            .mix(SUMMARIZE_COUNTS.out.read_tracking)
            .mix(SUMMARIZE_COUNTS.out.otu_summary)
    }

    COLLECT_RESULTS(
        ch_misc.collect(),
        ch_cluster.collect()
    )
}

// =============================================================================
// ENTRY POINT  —  help, run summary, USEARCH fetch, then CECILIA
// =============================================================================
workflow {

    if (Utils.asBool(params.help)) {
        log.info """
    ╔═══════════════════════════════════════════╗
    ║            C E C I L I A  v2.0            ║
    ║   Amplicon pipeline for Illumina data     ║
    ╚═══════════════════════════════════════════╝

    Usage:
        nextflow run cecilia.nf [options]
        sbatch run_cecilia.sb [options]

    ── Input / Output ───────────────────────────────────────────────────────────
      --rawdata          DIR     Input directory with *.fastq.gz files
                                 [default: ../rawdata]
      --outdir           DIR     Results directory
                                 [default: ../outputs]

    ── Tools & Databases ────────────────────────────────────────────────────────
      --usearch          PATH    USEARCH binary. Downloaded automatically on
                                 first run (v11.0.667, public domain CC0-1.0)
                                 [default: <pipeline>/.tools/usearch11.0.667_i86linux64]
      --usearch_fetch    BOOL    Allow that download; false = fail if absent
                                 [default: true]
      --sintax_db        PATH    FASTA database for SINTAX taxonomy (e.g. SILVA)
      --closedRef_db     PATH    FASTA database for closed-reference OTU picking

    ── DNA marker & read layout ─────────────────────────────────────────────────
      --dna_marker       STR     Marker label used in output names  [default: 16S]
      --paired           BOOL    true = paired-end R1+R2, false = single-end
                                 [default: true]

    ── Pipeline flow flags ───────────────────────────────────────────────────────
      --assemble         BOOL    Merge paired reads (USEARCH -fastq_mergepairs)
                                 [default: true]
      --rename           BOOL    Rename reads to sequential sample001/002/…
                                 [default: true]
      --count_start      INT     Number of the first renamed sample  [default: 01]
      --primers          BOOL    Strip primers with Cutadapt
                                 [default: true]
      --stripleft        BOOL    Trim fixed N bp from read 5' end
                                 [default: false]
      --stripleft_bp     INT     Bases to strip when stripleft=true  [default: 0]
      --test_length      BOOL    Run the 150–300 bp length sweep
                                 [default: false]
      --cluster_otu      BOOL    UPARSE OTU clustering               [default: true]
      --cluster_asv      BOOL    UNOISE3 ASV denoising               [default: true]
      --cluster_asv_to_otu BOOL  Cluster ASVs → OTUs at clustering_id[default: true]
      --cluster_swarm    BOOL    SWARM d=1 fastidious clustering      [default: true]
      --closed_ref_otu   BOOL    Closed-reference OTU picking         [default: true]
      --sintax_taxonomy  BOOL    SINTAX taxonomy assignment           [default: true]
      --summary_table    BOOL    Per-sample read tracking + OTU summary [default: true]

    ── Primer sequences ─────────────────────────────────────────────────────────
      --fwd_primer       STR     Forward primer sequence
                                 [default: GTGCCAGCMGCCGCGGTAA  515F Parada/Caporaso]
      --rev_primer       STR     Reverse primer sequence
                                 [default: GGACTACHVGGGTWTCTAAT 806R Apprill/Caporaso]
      --rev_primer_rc    STR     Reverse primer reverse complement
                                 [default: ATTAGAWACCCBDGTAGTCC]

    ── Quality filtering ────────────────────────────────────────────────────────
      --max_eerr         FLOAT   Max expected errors (USEARCH -fastq_filter)
                                 [default: 1.0]
      --filter_lengths   LIST    Truncation lengths, e.g. [200,225,250]
                                 [default: [200, 225, 250]]
      --user_length      INT     Extra truncation length             [default: 290]
      --min_overlap      INT     Min overlap for paired assembly; 0 = USEARCH
                                 default (16 bp)                     [default: 0]

    ── Clustering ───────────────────────────────────────────────────────────────
      --otu_min_size     INT     Min OTU size (UPARSE)               [default: 1]
      --zotu_min_size    INT     Min ZOTU size (UNOISE3)             [default: 8]
      --clustering_id    FLOAT   Identity for ASV→OTU clustering     [default: 0.97]

    ── Taxonomy ─────────────────────────────────────────────────────────────────
      --sintax_cutoff    FLOAT   SINTAX confidence cutoff            [default: 0.8]

    ── SLURM ────────────────────────────────────────────────────────────────────
      --slurm_account    STR     SLURM billing account               [default: glbrc]

    ─────────────────────────────────────────────────────────────────────────────
    Examples:
      # Fresh run with defaults
      sbatch run_cecilia.sb

      # ITS dataset with custom primers, skip SWARM
      sbatch run_cecilia.sb \\
        --dna_marker ITS \\
        --fwd_primer CTTGGTCATTTAGAGGAAGTAA \\
        --rev_primer TCCTCCGCTTATTGATATGC \\
        --rev_primer_rc GCATATCAATAAGCGGAGGA \\
        --cluster_swarm false \\
        --sintax_db /path/to/UNITE.fasta

      # Resume a previous run
      sbatch run_cecilia.sb -resume <session-id>
    ─────────────────────────────────────────────────────────────────────────────
        """.stripIndent()
        exit 0
    }

    // Unmerged pairs are not supported downstream: renaming would treat R1 and
    // R2 as two samples, and Cutadapt has no paired mode here. Stop early.
    if (Utils.asBool(params.paired) && !Utils.asBool(params.assemble)) {
        error """
        --paired true needs --assemble true: CECILIA does not process unmerged
        read pairs. If your pairs do not merge (amplicon longer than about
        2 x read length - 20 bp), run on R1 only: put only the R1 files in
        --rawdata and use --paired false.
        """.stripIndent()
    }

    // Put the USEARCH binary in place before any task is submitted. Downloads
    // it on first run; a no-op once .tools/ is populated. See
    // Utils.ensureUsearch() for why this cannot be a Nextflow process.
    def usearch_bin = Utils.ensureUsearch(
        params.usearch,
        params.usearch_url,
        params.usearch_sha256,
        Utils.asBool(params.usearch_fetch)
    )

    log.info """
    ╔═══════════════════════════════════════════╗
    ║            C E C I L I A  v2.0            ║
    ║   Amplicon pipeline for Illumina data     ║
    ╚═══════════════════════════════════════════╝

    usearch       : ${usearch_bin}
    rawdata       : ${params.rawdata}
    outdir        : ${params.outdir}
    dna_marker    : ${params.dna_marker}
    paired        : ${params.paired}
    assemble      : ${params.assemble}
    rename        : ${params.rename}
    primers       : ${params.primers}
    test_length   : ${params.test_length}
    cluster_otu   : ${params.cluster_otu}
    cluster_asv   : ${params.cluster_asv}
    cluster_swarm : ${params.cluster_swarm}
    closed_ref    : ${params.closed_ref_otu}
    sintax        : ${params.sintax_taxonomy}
    summary_table : ${params.summary_table}
    """.stripIndent()

    CECILIA()
}
