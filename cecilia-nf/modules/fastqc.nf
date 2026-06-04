// =============================================================================
// FASTQC  (step 02, per sample)
// Runs FastQC on each sample's decompressed read files.
// Also runs USEARCH -fastx_info per file (binary mounted from host).
// Parallel branch — does not block the main pipeline DAG.
// =============================================================================
process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/02_rawQuality_fqc",          mode: 'copy', pattern: "*.{html,zip}"
    publishDir "${params.outdir}/02_rawQuality_fqc/fastq_info", mode: 'copy', pattern: "*.info"
    publishDir "${params.outdir}/stats",                        mode: 'copy', pattern: "*.counts"

    input:
    tuple val(sample_id), path(reads)

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

    # Raw read count per file
    counts_file="${sample_id}.raw.counts"
    for f in ${reads.join(' ')}; do
        count=\$(( \$(wc -l < "\$f") / 4 ))
        echo "\$(basename \$f) : \$count" >> "\$counts_file"
    done
    """
}
