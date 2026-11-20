// =============================================================================
// CLOSED_REF_OTU  (step 14, optional: params.closed_ref_otu = true)
// USEARCH -closed_ref maps dereplicated reads to params.closedRef_db at
// params.clustering_id identity.  Useful for SynCom or non-overlapping
// amplicons; see USEARCH docs for known limitations.
// Input: tuple(len, uniques.fasta, pooled.fastq) from DEREPLICATE + POOL_READS
//
// Two stages, for a reason:
//
//   1. -closed_ref on the dereplicated uniques determines WHICH reference
//      sequences are detected at this truncation length.
//   2. -otutab on pooled.fastq determines HOW MANY reads each sample
//      contributes to each of them.
//
// Stage 2 is not optional bookkeeping. -fastx_uniques collapses identical
// sequences ACROSS samples and keeps a single representative label, so a
// unique reading `>sample001.36;size=16671;` carries reads from every sample
// but names only one. Building the sample x OTU table straight off the uniques
// therefore credits each sequence's entire abundance to one arbitrary sample —
// in practice the deepest, first-concatenated one. Re-mapping pooled.fastq,
// where every read still carries its own sampleNNN label, is what makes the
// per-sample columns real. This mirrors how CLUSTER_UPARSE / CLUSTER_UNOISE /
// CLUSTER_SWARM already build their tables.
// =============================================================================
process CLOSED_REF_OTU {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(uniques_fasta), path(pooled_fastq)

    output:
    tuple val(len), path("otutable_${len}bp_closedRef.txt"),   emit: table
    path("extra/otutable_${len}bp_closedRef.tabbed"),          emit: extra
    path("extra/closedRef_hits_${len}bp.fasta"),               emit: hits

    script:
    """
    set -euo pipefail
    mkdir -p extra

    # --- Stage 1: which references are present at this truncation length? ----
    ${params.usearch} \\
        -closed_ref ${uniques_fasta} \\
        -id ${params.clustering_id} \\
        -db ${params.closedRef_db} \\
        -strand both \\
        -otutabout extra/detected_refs_${len}bp.txt \\
        -tabbedout extra/otutable_${len}bp_closedRef.tabbed

    # OTU IDs in that table are the reference sequence labels that were hit.
    tail -n +2 extra/detected_refs_${len}bp.txt | cut -f1 | sed '/^\$/d' \\
        > extra/hit_ids_${len}bp.txt

    if [[ -s extra/hit_ids_${len}bp.txt ]]; then
        # --- Stage 2: per-sample counts from the labelled pooled reads ------
        ${params.usearch} \\
            -fastx_getseqs ${params.closedRef_db} \\
            -labels extra/hit_ids_${len}bp.txt \\
            -fastaout extra/closedRef_hits_${len}bp.fasta

        ${params.usearch} \\
            -threads ${task.cpus} \\
            -otutab ${pooled_fastq} \\
            -otus extra/closedRef_hits_${len}bp.fasta \\
            -id ${params.clustering_id} \\
            -strand both \\
            -otutabout otutable_${len}bp_closedRef.txt
    else
        echo "WARNING: no reference sequences detected at ${len} bp" >&2
        : > extra/closedRef_hits_${len}bp.fasta
        printf '#OTU ID\\n' > otutable_${len}bp_closedRef.txt
    fi
    """
}
