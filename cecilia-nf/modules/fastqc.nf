// =============================================================================
// FASTQC  (step 02)
// Runs FastQC and USEARCH fastx_info on all pooled samples for one read
// direction at a time (R1 forward reads together, R2 reverse reads together).
// Two SLURM jobs for paired data; one job for single-end.
// Parallel branch — does not block the main pipeline DAG.
// =============================================================================
process FASTQC {
    tag "${group}"
    publishDir "${params.outdir}/02_rawQuality_fqc",            mode: 'copy', pattern: "*.{html,zip}"
    publishDir "${params.outdir}/02_rawQuality_fqc/fastq_info", mode: 'copy', pattern: "*.info"
    publishDir "${params.outdir}/stats",                         mode: 'copy', pattern: "*.counts"

    input:
    tuple val(group), path(reads)   // 'R1', 'R2', or 'all' — all samples for that direction

    output:
    path("*.html"),   emit: html
    path("*.zip"),    emit: zip
    path("*.info"),   emit: info,   optional: true
    path("*.counts"), emit: counts

    script:
    """
    set -euo pipefail

    # Redirect JVM temp I/O to local /tmp to avoid SIGBUS on parallel filesystems
    export _JAVA_OPTIONS="-Djava.io.tmpdir=/tmp -Xss512k"

    fastqc --threads ${task.cpus} --dir /tmp --outdir . ${reads.join(' ')}

    for f in ${reads.join(' ')}; do
        base=\$(basename "\$f" .fastq)
        ${params.usearch} -fastx_info "\$f" -output "\${base}.info" 2>/dev/null || true
    done

    for f in ${reads.join(' ')}; do
        count=\$(( \$(wc -l < "\$f") / 4 ))
        echo "\$(basename \$f) : \$count"
    done > ${group}.raw.counts
    """
}
