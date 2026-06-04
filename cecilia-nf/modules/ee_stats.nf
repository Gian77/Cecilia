// =============================================================================
// EE_STATS  (step 07)
// Generates USEARCH expected-error statistics on the pooled read set.
// If params.stripleft = true:  first strips stripleft_bp 5' bases and emits
//                              trimmed.fastq for downstream steps.
// If params.stripleft = false: eestats2 runs directly on pooled.fastq.
// Also runs compareSeqQual.sh to flag reads where seq/quality lengths differ.
// =============================================================================
process EE_STATS {
    publishDir "${params.outdir}/07_readTrimDemux_usearch",    mode: 'copy', pattern: "*.{results,fastq}"
    publishDir "${params.outdir}/06_primerStrip_cutadapt",     mode: 'copy', pattern: "sequence_quality_nomatch.counts"
    publishDir "${params.outdir}/stats",                        mode: 'copy', pattern: "*.counts"

    input:
    path(pooled_fastq)

    output:
    path("eestats2.results"),                   emit: eestats
    path("trimmed.fastq"),                      emit: trimmed,  optional: true
    path("sequence_quality_nomatch.counts"),    emit: seqqual,  optional: true
    path("*.counts"),                           emit: counts,   optional: true

    script:
    if (params.stripleft) {
        """
        set -euo pipefail

        # Strip 5' bases
        ${params.usearch} \\
            -fastq_filter ${pooled_fastq} \\
            -fastq_stripleft ${params.stripleft_bp} \\
            -fastqout trimmed.fastq

        # Sequence/quality length sanity check on input directory
        compareSeqQual.sh .

        # EE statistics on trimmed reads
        ${params.usearch} \\
            -fastq_eestats2 trimmed.fastq \\
            -output eestats2.results \\
            -length_cutoffs 100,500,1

        # Per-sample read counts from trimmed headers (@sample001.N format)
        sed -n '1~4p' trimmed.fastq \\
            | sed 's/^@//' \\
            | awk -F'.' '{print \$1}' \\
            | sort | uniq -c \\
            | awk '{printf "%s : %s\\n", \$2, \$1}' > trimmed.counts
        """
    } else {
        """
        set -euo pipefail

        # Sequence/quality length sanity check
        compareSeqQual.sh .

        # EE statistics on pooled reads
        ${params.usearch} \\
            -fastq_eestats2 ${pooled_fastq} \\
            -output eestats2.results \\
            -length_cutoffs 100,500,1

        # Per-sample read counts from pooled headers (@sample001.N format)
        sed -n '1~4p' ${pooled_fastq} \\
            | sed 's/^@//' \\
            | awk -F'.' '{print \$1}' \\
            | sort | uniq -c \\
            | awk '{printf "%s : %s\\n", \$2, \$1}' > pooled.counts
        """
    }
}
