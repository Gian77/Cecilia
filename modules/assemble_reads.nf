// =============================================================================
// ASSEMBLE_READS  (step 04, optional: params.assemble = true, per sample)
// USEARCH -fastq_mergepairs merges paired-end reads into one FASTQ per sample.
// =============================================================================
process ASSEMBLE_READS {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_readAssembly_usearch", mode: 'copy', pattern: "*.fastq"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*assembled.fastq"), emit: reads
    path("*.counts"),                               emit: counts

    script:
    def r1  = reads[0]
    def r2  = reads[1]
    def out = "${sample_id}_assembled.fastq"
    // --min_overlap 0 keeps USEARCH's own default (16 bp)
    def min_ovl = params.min_overlap.toString().toInteger()
    def ovl_flag = min_ovl > 0 ? "-fastq_minovlen ${min_ovl}" : ""
    """
    set -euo pipefail

    ${params.usearch} \\
        -fastq_mergepairs ${r1} \\
        -reverse ${r2} \\
        ${ovl_flag} \\
        -fastqout ${out}

    count=\$(( \$(wc -l < "${out}") / 4 ))
    echo "${out} : \$count" > "${sample_id}.assembled.counts"
    """
}
