// =============================================================================
// FASTQC  (step 02, all samples together)
// Runs FastQC on all decompressed read files in a single job.
// Also runs USEARCH -fastx_info per file (binary mounted from host).
// Parallel branch — does not block the main pipeline DAG.
// =============================================================================
process FASTQC {
    publishDir "${params.outdir}/02_rawQuality_fqc",            mode: 'copy', pattern: "*.{html,zip}"
    publishDir "${params.outdir}/02_rawQuality_fqc/fastq_info", mode: 'copy', pattern: "*.info"
    publishDir "${params.outdir}/stats",                         mode: 'copy', pattern: "*.counts"

    input:
    path(reads)   // all decompressed read files collected into one invocation

    output:
    path("*.html"),   emit: html
    path("*.zip"),    emit: zip
    path("*.info"),   emit: info,   optional: true
    path("*.counts"), emit: counts

    script:
    """
    set -euo pipefail

    # Redirect JVM temp I/O to local /tmp to avoid SIGBUS on parallel filesystems
    # (GPFS/BeeGFS mmap behaviour triggers bus errors in the OpenJDK G1GC).
    export _JAVA_OPTIONS="-Djava.io.tmpdir=/tmp -Xss512k"

    # FastQC: --dir /tmp keeps all temp files off the parallel filesystem
    fastqc --threads ${task.cpus} --dir /tmp --outdir . ${reads.join(' ')}

    # USEARCH fastx_info per file (USEARCH binary is bind-mounted from host)
    for f in ${reads.join(' ')}; do
        base=\$(basename "\$f" .fastq)
        ${params.usearch} -fastx_info "\$f" -output "\${base}.info" 2>/dev/null || true
    done

    # Raw read counts — one line per file, all in a single summary file
    for f in ${reads.join(' ')}; do
        count=\$(( \$(wc -l < "\$f") / 4 ))
        echo "\$(basename \$f) : \$count"
    done > raw.counts
    """
}
