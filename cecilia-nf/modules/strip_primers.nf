// =============================================================================
// STRIP_PRIMERS  (step 06, per sample)
// Cutadapt removes primers/adapters. Adapter flags depend on assembly state:
//   assembled  → -g fwd_primer -a rev_primer_rc  -n 2   (both orientations)
//   unassembled, R1+R2 → -g fwd_primer on R1 (common single-file case post-rename)
//   unassembled, R1-only → -g fwd_primer
// If params.primers = false, the file is passed through unchanged.
// =============================================================================
process STRIP_PRIMERS {
    tag "${sample_id}"
    publishDir "${params.outdir}/06_primerStrip_cutadapt", mode: 'copy', pattern: "*.fastq"

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*stripped.fastq"), emit: fastq
    path("*.counts"),                              emit: counts

    script:
    def out = "${sample_id}_stripped.fastq"
    if (params.primers) {
        def adapter_flags = params.assemble
            ? "-g ${params.fwd_primer} -a ${params.rev_primer_rc} -n 2"
            : "-g ${params.fwd_primer}"
        """
        set -euo pipefail

        cutadapt ${adapter_flags} \\
            -j ${task.cpus} \\
            -e 0.01 \\
            --discard-untrimmed \\
            --match-read-wildcards \\
            -o ${out} \\
            ${fastq}

        count=\$(( \$(wc -l < "${out}") / 4 ))
        echo "${out} : \$count" > "${sample_id}.stripped.counts"
        """
    } else {
        """
        set -euo pipefail
        cp ${fastq} ${out}
        count=\$(( \$(wc -l < "${out}") / 4 ))
        echo "${out} : \$count" > "${sample_id}.stripped.counts"
        """
    }
}
