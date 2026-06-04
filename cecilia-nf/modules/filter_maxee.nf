// =============================================================================
// FILTER_MAXEE  (step 09, one process per truncation length — all parallel)
// Filters pooled/trimmed reads at params.max_eerr and the given bp length.
// Outputs nothing and exits cleanly if no reads survive the filter.
// Runs FastQC on the filtered file for post-filter quality review.
//
// Input  : tuple(len, source_fastq)  — built in workflow with Channel.combine()
// Output : tuple(len, filtered.fastq) — declared optional; suppressed if empty
// =============================================================================
process FILTER_MAXEE {
    tag "len=${len}"
    publishDir "${params.outdir}/09_filterMaxee_usearch", mode: 'copy'

    input:
    tuple val(len), path(source_fastq)

    output:
    tuple val(len), path("filtered_${len}bp.fastq"), emit: filtered,     optional: true
    path("*.html"),                                   emit: fastqc_html,  optional: true
    path("*.zip"),                                    emit: fastqc_zip,   optional: true

    script:
    """
    set -euo pipefail

    ${params.usearch} \\
        -fastq_filter ${source_fastq} \\
        -fastq_maxee ${params.max_eerr} \\
        -fastq_trunclen ${len} \\
        -fastq_maxns 0 \\
        -fastqout filtered_${len}bp.fastq

    if [[ ! -s filtered_${len}bp.fastq ]]; then
        echo "WARNING: no reads survived filter at ${len} bp — skipping FastQC" >&2
        rm -f filtered_${len}bp.fastq
        exit 0
    fi

    export _JAVA_OPTIONS="-Djava.io.tmpdir=/tmp -Xss512k"
    fastqc filtered_${len}bp.fastq \\
        --threads ${task.cpus} \\
        --dir /tmp \\
        --outdir .
    """
}
