// =============================================================================
// RENAME_READS  (step 05, optional: params.rename = true)
// Collects ALL per-sample FASTQs into one job (matching the original behaviour
// of renameSeq.sh, which assigns sequential numbers across the full dataset).
//
// Input  : all per-sample FASTQs, collected via .map { id,f -> f }.collect()
// Output : Sample*.fastq files + file_mapping.txt
//          Downstream channels re-parallelise with .flatten().map{ f->[f.baseName,f] }
//
// Note: uses inline awk rather than calling bin/renameSeq.sh so that cp (not mv)
// is used — Nextflow stages inputs as symlinks, which cannot be renamed in-place.
// =============================================================================
process RENAME_READS {
    publishDir "${params.outdir}/05_readRenamed_bash", mode: 'copy'

    input:
    path(fastqs)   // collected list from all upstream samples

    output:
    path("sample*.fastq"), emit: reads
    path("file_mapping.txt"), emit: mapping

    script:
    """
    set -euo pipefail

    count=${params.count_start as Integer}

    # Sort for deterministic sequential numbering
    for f in \$(ls *.fastq | sort); do
        [ -f "\$f" ] || continue
        count_padded=\$(printf "%03d" "\$count")
        sample_name="sample\${count_padded}"
        new_name="\${sample_name}.fastq"

        # Record mapping (original → new)
        echo -e "\$f\t\$new_name" >> file_mapping.txt

        # Copy (not mv) because NF stages inputs as symlinks
        cp "\$f" "\$new_name"

        # Rewrite FASTQ @header lines: @sample001.1, @sample001.2, ...
        awk -v prefix="\$sample_name" '
            { if (NR % 4 == 1) {
                read_count++
                printf "@%s.%d\\n", prefix, read_count
              } else {
                print
              } }' "\$new_name" > tmp_rename.fastq
        mv tmp_rename.fastq "\$new_name"

        count=\$(( count + 1 ))
    done
    """
}
