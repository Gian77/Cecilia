// =============================================================================
// CLUSTER_UNOISE  (step 12, optional: params.cluster_asv = true)
// USEARCH -unoise3 → ASV FASTA (Zotu headers rewritten to asv_NNN)
//                 + ASV abundance table mapped against pooled.fastq
// Input: tuple(len, uniques.fasta, pooled.fastq) — combined in workflow
// =============================================================================
process CLUSTER_UNOISE {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(uniques_fasta), path(pooled_fastq)

    output:
    tuple val(len), path("asv_${len}bp.fasta"),        emit: asv
    path("asvtable_UNOISE_${len}bp.txt"),              emit: table
    path("extra/unoise_asv_${len}bp.txt"),             emit: extra

    script:
    """
    set -euo pipefail
    mkdir -p extra

    ${params.usearch} \\
        -unoise3 ${uniques_fasta} \\
        -minsize ${params.zotu_min_size} \\
        -tabbedout extra/unoise_asv_${len}bp.txt \\
        -zotus asv_${len}bp.fasta

    # Rewrite >ZotuN headers to >asv_NNN
    awk '/^>/ {printf(">asv_%03d\\n", ++i); next} {print}' \\
        asv_${len}bp.fasta > asv_${len}bp.tmp
    mv asv_${len}bp.tmp asv_${len}bp.fasta

    # usearch -otutab segfaults (exit 139) on an empty -zotus FASTA, so write a
    # header-only table instead, as CLOSED_REF_OTU does.
    if [[ -s asv_${len}bp.fasta ]]; then
        ${params.usearch} \\
            -threads ${task.cpus} \\
            -otutab ${pooled_fastq} \\
            -zotus asv_${len}bp.fasta \\
            -otutabout asvtable_UNOISE_${len}bp.txt
    else
        echo "WARNING: UNOISE3 produced no ASVs at ${len} bp" >&2
        printf '#OTU ID\\n' > asvtable_UNOISE_${len}bp.txt
    fi
    """
}
