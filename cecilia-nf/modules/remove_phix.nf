// =============================================================================
// REMOVE_PHIX  (step 03, per sample)
// USEARCH -filter_phix on each sample.
// Paired-end (params.paired = true AND params.assemble = true):
//   two output files: *_R1_nophix.fastq, *_R2_nophix.fastq
// Single-end R1-only or post-assembly single file:
//   one output file: *_nophix.fastq
// =============================================================================
process REMOVE_PHIX {
    tag "${sample_id}"
    publishDir "${params.outdir}/03_removePhix_usearch", mode: 'copy', pattern: "*.fastq"
    publishDir "${params.outdir}/stats",                 mode: 'copy', pattern: "*.counts"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*nophix.fastq"), emit: reads
    path("*.counts"),                            emit: counts

    script:
    def r1   = reads[0]
    def r2   = reads.size() > 1 ? reads[1] : null
    def paired_cmd = r2
        ? """
          ${params.usearch} \\
              -filter_phix ${r1} \\
              -reverse ${r2} \\
              -threads ${task.cpus} \\
              -output ${r1.baseName}_nophix.fastq \\
              -output2 ${r2.baseName}_nophix.fastq
          """
        : """
          ${params.usearch} \\
              -filter_phix ${r1} \\
              -threads ${task.cpus} \\
              -output ${r1.baseName}_nophix.fastq
          """
    """
    set -euo pipefail

    ${paired_cmd}

    # Read counts after PhiX removal
    counts_file="${sample_id}.nophix.counts"
    for f in *nophix.fastq; do
        count=\$(( \$(wc -l < "\$f") / 4 ))
        echo "\$f : \$count" >> "\$counts_file"
    done
    """
}
