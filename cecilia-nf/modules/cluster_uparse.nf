// =============================================================================
// CLUSTER_UPARSE  (step 11, optional: params.cluster_otu = true)
// USEARCH -cluster_otus → OTU FASTA (relabelled otu_NNN)
//                       + OTU abundance table mapped against pooled.fastq
// Input: tuple(len, uniques.fasta, pooled.fastq) — combined in workflow
// =============================================================================
process CLUSTER_UPARSE {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(uniques_fasta), path(pooled_fastq)

    output:
    tuple val(len), path("otus_${len}bp.fasta"),       emit: fasta
    path("otutable_UPARSE_${len}bp.txt"),              emit: table
    path("extra/uparse_otus_${len}bp.txt"),            emit: extra

    script:
    """
    set -euo pipefail
    mkdir -p extra

    ${params.usearch} \\
        -minsize ${params.otu_min_size} \\
        -relabel otu_ \\
        -cluster_otus ${uniques_fasta} \\
        -otus otus_${len}bp.fasta \\
        -uparseout extra/uparse_otus_${len}bp.txt

    ${params.usearch} \\
        -threads ${task.cpus} \\
        -otutab ${pooled_fastq} \\
        -otus otus_${len}bp.fasta \\
        -otutabout otutable_UPARSE_${len}bp.txt
    """
}
