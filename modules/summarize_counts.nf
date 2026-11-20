// =============================================================================
// SUMMARIZE_COUNTS  (step 19, optional: params.summary_table = true)
// Joins the per-stage counts the pipeline already emits into two tidy TSVs:
//   read_tracking.tsv  one row per sample  (raw, nophix, assembled, stripped)
//   otu_summary.tsv    long format, one row per (sample, length, method)
//
// Deliberately self-contained: it reads the raw FASTQs itself rather than
// having DECOMPRESS emit raw counts. Adding an output to DECOMPRESS would
// change that process's script hash and invalidate every downstream task in
// the -resume cache, costing a full rerun on existing analyses.
//
// All optional inputs arrive through the single collected `aux_files` channel,
// so absent stages (no rename, no assembly) simply stage no file and the
// script falls back — no placeholder plumbing needed.
// =============================================================================
process SUMMARIZE_COUNTS {
    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    path(rawdir)
    path(aux_files)

    output:
    path("read_tracking.tsv"), emit: read_tracking
    path("otu_summary.tsv"),   emit: otu_summary

    script:
    """
    set -euo pipefail
    summarizeCounts.sh ${rawdir} .
    """
}
