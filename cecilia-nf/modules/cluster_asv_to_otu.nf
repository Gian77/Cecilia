// =============================================================================
// CLUSTER_ASV_TO_OTU  (step 13, optional: params.cluster_asv_to_otu = true)
// Clusters ASV FASTA at params.clustering_id (default 97%) using
// USEARCH -sortbylength → -cluster_smallmem → OTU table.
// Requires params.cluster_asv = true (runs on CLUSTER_UNOISE output).
// Input: tuple(len, asv.fasta, pooled.fastq) — combined in workflow
// =============================================================================
process CLUSTER_ASV_TO_OTU {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(asv_fasta), path(pooled_fastq)

    output:
    tuple val(len), path("asv_${len}bp_to97otus.fasta"),  emit: fasta
    path("otutable_asv_${len}bp_to97otus.txt"),           emit: table
    path("extra/asv_${len}bp_to97otus.clust"),            emit: extra

    script:
    """
    set -euo pipefail
    mkdir -p extra

    ${params.usearch} \\
        -sortbylength ${asv_fasta} \\
        -fastaout extra/asv_${len}bp_sorted.fasta

    ${params.usearch} \\
        -cluster_smallmem extra/asv_${len}bp_sorted.fasta \\
        -id ${params.clustering_id} \\
        -centroids asv_${len}bp_to97otus.fasta \\
        -uc extra/asv_${len}bp_to97otus.clust

    ${params.usearch} \\
        -threads ${task.cpus} \\
        -otutab ${pooled_fastq} \\
        -otus asv_${len}bp_to97otus.fasta \\
        -otutabout otutable_asv_${len}bp_to97otus.txt
    """
}
