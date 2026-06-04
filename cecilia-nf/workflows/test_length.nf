nextflow.enable.dsl = 2

include { TEST_LENGTH } from '../modules/test_length'

// =============================================================================
// Optional subworkflow: sweep read lengths 150–300 bp and plot OTU yield
// Triggered only when params.test_length = true
// =============================================================================
workflow TEST_LENGTH_WF {

    take:
    ch_source_fastq   // pooled.fastq or trimmed.fastq

    main:
    TEST_LENGTH(ch_source_fastq)

    emit:
    results = TEST_LENGTH.out.results
    plots   = TEST_LENGTH.out.plots
}
