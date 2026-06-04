// =============================================================================
// DEREPLICATE  (step 10, one process per filtered file — all parallel)
// USEARCH -fastx_uniques → uniques_${len}bp.fasta  (UPARSE / UNOISE / closed-ref)
// SWARM 0-diff pre-clustering → uniques_${len}bp_linear.fasta  (SWARM only)
// The linear output is optional: only built when params.cluster_swarm = true.
// =============================================================================
process DEREPLICATE {
    tag "len=${len}"
    publishDir "${params.outdir}/10_dereplicateReads_usearch", mode: 'copy'

    input:
    tuple val(len), path(filtered_fastq)

    output:
    tuple val(len), path("uniques_${len}bp.fasta"),        emit: uniques
    tuple val(len), path("uniques_${len}bp_linear.fasta"), emit: linear, optional: true

    script:
    def swarm_prep = params.cluster_swarm ? """
        # Linearise FASTA for SWARM (awk one-liner from original pipeline)
        awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0); N++; next} \\
             {printf("%s",\$0)} END {printf("\\n")}' \\
            uniques_${len}bp.fasta \\
            | tr "\\t" "\\n" \\
            | sed -e 's/\\( \\).*\\(;.\\)/\\1\\2/' \\
            | sed 's/ //' \\
            | sed 's/.\$//' \\
            > uniques_${len}bp_linear.temp

        swarm \\
            --threads ${task.cpus} \\
            --differences 0 \\
            -w uniques_${len}bp_linear.fasta \\
            -z \\
            -o /dev/null \\
            uniques_${len}bp_linear.temp

        rm -f uniques_${len}bp_linear.temp
    """ : ""
    """
    set -euo pipefail

    ${params.usearch} \\
        -threads ${task.cpus} \\
        -fastx_uniques ${filtered_fastq} \\
        -fastaout uniques_${len}bp.fasta \\
        -sizeout

    ${swarm_prep}
    """
}
