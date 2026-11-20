// =============================================================================
// BUILD_SINTAX_DB  (step 17a)
// Builds a USEARCH .udb index from the SINTAX reference FASTA once.
// Decompresses with gunzip -k (keep original) if the file is gzip-compressed.
// The .udb is shared across all TAXONOMY_SINTAX process calls via .combine().
// =============================================================================
process BUILD_SINTAX_DB {
    // Deliberately NOT published. The .udb is a derived index, not a result:
    // for SILVA 138 it is ~13 GB, and `mode: 'copy'` duplicated it into every
    // run's outdir (13 of the 15 GB in a test run's outputs/). Nextflow keeps
    // it in work/ for downstream tasks and for -resume, and it rebuilds from
    // params.sintax_db in minutes. Publishing it only burns $HOME quota.

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
