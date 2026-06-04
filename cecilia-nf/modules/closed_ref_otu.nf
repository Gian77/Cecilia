// =============================================================================
// CLOSED_REF_OTU  (step 14, optional: params.closed_ref_otu = true)
// USEARCH -closed_ref maps dereplicated reads to params.closedRef_db at
// params.clustering_id identity.  Useful for SynCom or non-overlapping
// amplicons; see USEARCH docs for known limitations.
// Input: tuple(len, uniques.fasta) from DEREPLICATE
// =============================================================================
process CLOSED_REF_OTU {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(uniques_fasta)

    output:
    tuple val(len), path("otutable_${len}bp_closedRef.txt"),   emit: table
    path("extra/otutable_${len}bp_closedRef.tabbed"),          emit: extra

    script:
    """
    set -euo pipefail
    mkdir -p extra

    ${params.usearch} \\
        -closed_ref ${uniques_fasta} \\
        -id ${params.clustering_id} \\
        -db ${params.closedRef_db} \\
        -strand both \\
        -otutabout otutable_${len}bp_closedRef.txt \\
        -tabbedout extra/otutable_${len}bp_closedRef.tabbed
    """
}
