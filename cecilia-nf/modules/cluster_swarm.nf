// =============================================================================
// CLUSTER_SWARM  (step 15, optional: params.cluster_swarm = true)
// SWARM --differences 1 --fastidious on the 0-diff pre-clustered FASTA
// from DEREPLICATE, then USEARCH -unoise3 on the SWARM seeds to produce
// the final SWARM FASTA (headers rewritten to swarm_NNN) + OTU table.
// Input: tuple(len, linear.fasta, pooled.fastq) — combined in workflow
// =============================================================================
process CLUSTER_SWARM {
    tag "len=${len}"
    publishDir "${params.outdir}/11_clustered_otu_asv_usearch", mode: 'copy'

    input:
    tuple val(len), path(linear_fasta), path(pooled_fastq)

    output:
    tuple val(len), path("swarms_${len}bp.fasta"),       emit: fasta
    path("otutable_SWARM_${len}bp.txt"),                 emit: table
    path("extra/clusters_${len}bp_stats.txt"),           emit: stats
    path("extra/clusters_${len}bp.fasta"),               emit: seeds

    script:
    """
    set -euo pipefail
    mkdir -p extra

    swarm \\
        --differences 1 \\
        --fastidious \\
        --usearch-abundance \\
        --threads ${task.cpus} \\
        --statistics-file extra/clusters_${len}bp_stats.txt \\
        --seeds extra/clusters_${len}bp.fasta \\
        ${linear_fasta} > /dev/null

    ${params.usearch} \\
        -threads ${task.cpus} \\
        -unoise3 extra/clusters_${len}bp.fasta \\
        -tabbedout extra/swarms_${len}bp.txt \\
        -zotus swarms_${len}bp.fasta

    # Rewrite >ZotuN headers to >swarm_NNN
    awk '/^>/ {printf(">swarm_%03d\\n", ++i); next} {print}' \\
        swarms_${len}bp.fasta > swarms_${len}bp.tmp
    mv swarms_${len}bp.tmp swarms_${len}bp.fasta

    ${params.usearch} \\
        -threads ${task.cpus} \\
        -otutab ${pooled_fastq} \\
        -zotus swarms_${len}bp.fasta \\
        -otutabout otutable_SWARM_${len}bp.txt
    """
}
