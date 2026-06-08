// =============================================================================
// TEST_LENGTH  (step 08, optional subworkflow — params.test_length = true)
// Sweeps truncation lengths 150–300 bp: filter → derep → cluster_otus at
// each length, records counts in testLength.results, then plots with R.
// Requires: USEARCH (bind-mounted), Rscript + tidyverse (from rocker/tidyverse).
// bin/testLength.R is on PATH via Nextflow's bin/ convention.
// =============================================================================
process TEST_LENGTH {
    publishDir "${params.outdir}/08_testLength_usearch", mode: 'copy'

    input:
    path(source_fastq)

    output:
    path("testLength.results"), emit: results
    path("*.pdf"),              emit: plots

    script:
    def stripleft_flag = params.stripleft ? "-fastq_stripleft ${params.stripleft_bp}" : ""
    """
    set -euo pipefail

    # Rename read headers to @newstring.N for consistent grep-counting
    awk 'NR%4==1{sub(/.*\\./, "@newstring.")} 1' ${source_fastq} > source_renamed.fastq

    echo -e "Length\tFiltered\tUniques\tOtus" > testLength.results

    for i in \$(seq -w 150 300); do

        ${params.usearch} -threads ${task.cpus} \\
            -fastq_filter source_renamed.fastq \\
            -fastq_maxee ${params.max_eerr} \\
            ${stripleft_flag} \\
            -fastq_trunclen \$i \\
            -fastq_maxns 0 \\
            -fastqout filtered_\${i}.fastq

        if [[ ! -s filtered_\${i}.fastq ]]; then
            echo -e "\${i}\t0\t0\t0" >> testLength.results
            rm -f filtered_\${i}.fastq
            continue
        fi

        ${params.usearch} -threads ${task.cpus} \\
            -fastx_uniques filtered_\${i}.fastq \\
            -fastaout uniques_\${i}.fasta \\
            -sizeout

        ${params.usearch} -threads ${task.cpus} \\
            -minsize 2 \\
            -relabel otu_ \\
            -cluster_otus uniques_\${i}.fasta \\
            -otus otus_\${i}.fasta \\
            -uparseout uparse_\${i}.txt 2>/dev/null || true

        filt=\$(grep -c "^@newstring"  filtered_\${i}.fastq  2>/dev/null || echo 0)
        uniq=\$(grep -c "^>newstring"  uniques_\${i}.fasta   2>/dev/null || echo 0)
        otu=\$(grep  -c "^>otu"        otus_\${i}.fasta      2>/dev/null || echo 0)

        echo -e "\${i}\t\${filt}\t\${uniq}\t\${otu}" >> testLength.results

        rm -f filtered_\${i}.fastq uniques_\${i}.fasta otus_\${i}.fasta uparse_\${i}.txt
    done

    rm -f source_renamed.fastq

    # Plot: testLength.R reads testLength.results from the working directory
    testLength.R
    """
}
