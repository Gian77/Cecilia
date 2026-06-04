// =============================================================================
// TAXONOMY_SINTAX  (step 17, optional: params.sintax_taxonomy = true)
// USEARCH -sintax on each cluster FASTA against the pre-built .udb index.
// One process per FASTA file — all run in parallel.
// Input: tuple(len, cluster.fasta, sintax.udb)
//   built in workflow with: ch_cluster_fastas.combine(BUILD_SINTAX_DB.out.udb)
// =============================================================================
process TAXONOMY_SINTAX {
    tag "len=${len}"
    publishDir "${params.outdir}/13_taxonomy_sintax", mode: 'copy'

    input:
    tuple val(len), path(cluster_fasta), path(sintax_udb)

    output:
    tuple val(len), path("${cluster_fasta.baseName}.sintax"), emit: sintax

    script:
    """
    set -euo pipefail

    ${params.usearch} \\
        -sintax ${cluster_fasta} \\
        -db ${sintax_udb} \\
        -tabbedout ${cluster_fasta.baseName}.sintax \\
        -strand both \\
        -sintax_cutoff ${params.sintax_cutoff} \\
        -threads ${task.cpus}
    """
}
