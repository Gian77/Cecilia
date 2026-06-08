#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// =============================================================================
// Cecilia-nf  —  main.nf
// usearCh basEd ampliCon pIpeLine for Illumina dAta
// Gian M.N. Benucci, Ph.D.  |  benucci@msu.edu
// =============================================================================

include { CECILIA } from './workflows/cecilia'

workflow {

    if (params.help) {
        log.info """
    ╔═══════════════════════════════════════════╗
    ║            C E C I L I A  v2.0            ║
    ║   Amplicon pipeline for Illumina data     ║
    ╚═══════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options]
        sbatch run_cecilia.sb [options]

    ── Input / Output ───────────────────────────────────────────────────────────
      --rawdata          DIR     Input directory with *.fastq.gz files
                                 [default: ../rawdata]
      --outdir           DIR     Results directory
                                 [default: ../outputs]

    ── Tools & Databases (site-specific paths) ──────────────────────────────────
      --usearch          PATH    USEARCH binary (bind-mounted into containers)
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
      --primers          BOOL    Strip primers with Cutadapt
                                 [default: true]
      --stripleft        BOOL    Trim fixed N bp from read 5' end
                                 [default: false]
      --stripleft_bp     INT     Bases to strip when stripleft=true  [default: 0]
      --test_length      BOOL    Run 150–300 bp length-sweep subworkflow
                                 [default: false]
      --cluster_otu      BOOL    UPARSE OTU clustering               [default: true]
      --cluster_asv      BOOL    UNOISE3 ASV denoising               [default: true]
      --cluster_asv_to_otu BOOL  Cluster ASVs → OTUs at clustering_id[default: true]
      --cluster_swarm    BOOL    SWARM d=1 fastidious clustering      [default: true]
      --closed_ref_otu   BOOL    Closed-reference OTU picking         [default: true]
      --sintax_taxonomy  BOOL    SINTAX taxonomy assignment           [default: true]

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
      --min_overlap      INT     Min overlap for paired assembly     [default: 0]

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

    log.info """
    ╔═══════════════════════════════════════════╗
    ║            C E C I L I A  v2.0            ║
    ║   Amplicon pipeline for Illumina data     ║
    ╚═══════════════════════════════════════════╝

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
    """.stripIndent()

    CECILIA()
}
