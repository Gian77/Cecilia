// =============================================================================
// POOL_READS  (step 06b)
// Collects all per-sample stripped FASTQs into a single pooled.fastq.
// Also creates a 500-read FASTA subset for QC, and counts reads per sample.
// This is a single job; inputs arrive via .collect() in the workflow.
// =============================================================================
process POOL_READS {
    publishDir "${params.outdir}/06_primerStrip_cutadapt", mode: 'copy'
    publishDir "${params.outdir}/stats",                   mode: 'copy', pattern: "pooled.counts"

    input:
    path(stripped_fastqs)   // collected list of *_stripped.fastq files

    output:
    path("pooled.fastq"),     emit: pooled_fastq
    path("subset_500.fasta"), emit: subset_fasta
    path("pooled.counts"),    emit: counts

    script:
    """
    set -euo pipefail

    # Pool all stripped reads
    cat *_stripped.fastq > pooled.fastq

    # 500-read FASTA subset for downstream QC
    seqtk sample -s100 pooled.fastq 500 > subset_500.fastq
    seqtk seq -aQ64 subset_500.fastq   > subset_500.fasta
    rm -f subset_500.fastq

    # Per-sample read counts extracted from pooled.fastq headers
    # Headers look like: @sample001.1 (written by RENAME_READS awk step)
    awk 'NR%4==1' pooled.fastq \\
        | sed 's/^@//' \\
        | awk -F'.' '{print \$1}' \\
        | sort | uniq -c \\
        | awk '{printf "%s : %s\\n", \$2, \$1}' > pooled.counts
    """
}
