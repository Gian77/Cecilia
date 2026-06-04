// =============================================================================
// COLLECT_RESULTS  (step 18)
// Gathers all tracked output files passed in from upstream channels,
// bundles them into a Cecilia-results/ directory, and creates a dated
// .tar.gz archive named after the read layout (assembled / paired / single).
// All publishDir writes are already done upstream; this step creates the
// portable archive for handoff.
// =============================================================================
process COLLECT_RESULTS {
    publishDir "${params.outdir}/14_getResults_bash", mode: 'copy'

    input:
    path(misc_files)     // QC, EE stats, subset FASTA, mapping, test-length outputs
    path(cluster_files)  // cluster FASTAs, OTU tables, sintax outputs

    output:
    path("Cecilia-results_*.tar.gz"), emit: archive

    script:
    def suffix = params.assemble ? 'assembled' : (params.paired ? 'paired' : 'single')
    """
    set -euo pipefail

    results_dir="Cecilia-results"
    mkdir -p "\${results_dir}"

    # Copy all tracked files — skip directories passed via collect()
    for f in ${misc_files.join(' ')} ${cluster_files.join(' ')}; do
        [ -f "\$f" ] && cp "\$f" "\${results_dir}/"
    done

    # Timestamp and package
    timestamp=\$(date +"%d-%m-%Y")
    archive="Cecilia-results_${suffix}_\${timestamp}.tar.gz"

    tar -zcf "\${archive}" "\${results_dir}/"

    echo "Archive created: \${archive}"
    ls -lh "\${archive}"
    """
}
