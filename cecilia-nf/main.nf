#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// =============================================================================
// Cecilia-nf  —  main.nf
// usearCh basEd ampliCon pIpeLine for Illumina dAta
// Gian M.N. Benucci, Ph.D.  |  benucci@msu.edu
// =============================================================================

include { CECILIA } from './workflows/cecilia'

workflow {

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
