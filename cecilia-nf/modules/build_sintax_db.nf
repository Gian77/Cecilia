// =============================================================================
// BUILD_SINTAX_DB  (step 17a)
// Builds a USEARCH .udb index from the SINTAX reference FASTA once.
// Decompresses with gunzip -k (keep original) if the file is gzip-compressed.
// The .udb is shared across all TAXONOMY_SINTAX process calls via .combine().
// =============================================================================
process BUILD_SINTAX_DB {
    publishDir "${params.outdir}/13_taxonomy_sintax", mode: 'copy'

    input:
    path(sintax_fasta)

    output:
    path("*.udb"), emit: udb

    script:
    def is_gz   = sintax_fasta.name.endsWith('.gz')
    def db_file = is_gz ? sintax_fasta.name.replace('.gz', '') : sintax_fasta.name
    def db_name = db_file.replaceAll(/\.fasta$/, '')
    """
    set -euo pipefail

    ${is_gz ? "gunzip -k ${sintax_fasta}" : ""}

    ${params.usearch} \\
        -makeudb_usearch ${db_file} \\
        -output ${db_name}.udb
    """
}
