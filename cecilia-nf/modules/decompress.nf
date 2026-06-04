// =============================================================================
// DECOMPRESS  (step 01, per sample)
// Decompresses .fastq.gz / .fastq.bz2 → .fastq; plain .fastq files are copied.
// Renames non-standard 1_1 / 1_2 suffixes to R1 / R2.
// =============================================================================
process DECOMPRESS {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_decompressed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fastq"), emit: reads

    script:
    """
    set -euo pipefail

    for f in ${reads.join(' ')}; do
        base=\$(basename "\$f")
        if [[ "\$f" == *.fastq.bz2 ]]; then
            bzip2 -cd "\$f" > "\${base%.bz2}"
        elif [[ "\$f" == *.fastq.gz ]]; then
            gzip -cd  "\$f" > "\${base%.gz}"
        elif [[ "\$f" == *.fastq ]]; then
            cp "\$f" "\$base"
        else
            echo "WARNING: \$f has unrecognised format — skipping" >&2
        fi
    done

    # Rename non-standard Illumina suffixes
    shopt -s nullglob
    for f in *1_1.fastq; do mv "\$f" "\${f//1_1.fastq/R1.fastq}"; done
    for f in *1_2.fastq; do mv "\$f" "\${f//1_2.fastq/R2.fastq}"; done
    """
}
