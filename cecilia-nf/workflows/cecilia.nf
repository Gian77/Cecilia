nextflow.enable.dsl = 2

include { DECOMPRESS        } from '../modules/decompress'
include { FASTQC            } from '../modules/fastqc'
include { REMOVE_PHIX       } from '../modules/remove_phix'
include { ASSEMBLE_READS    } from '../modules/assemble_reads'
include { RENAME_READS      } from '../modules/rename_reads'
include { STRIP_PRIMERS     } from '../modules/strip_primers'
include { POOL_READS        } from '../modules/pool_reads'
include { EE_STATS          } from '../modules/ee_stats'
include { FILTER_MAXEE      } from '../modules/filter_maxee'
include { DEREPLICATE       } from '../modules/dereplicate'
include { CLUSTER_UPARSE    } from '../modules/cluster_uparse'
include { CLUSTER_UNOISE    } from '../modules/cluster_unoise'
include { CLUSTER_ASV_TO_OTU } from '../modules/cluster_asv_to_otu'
include { CLOSED_REF_OTU    } from '../modules/closed_ref_otu'
include { CLUSTER_SWARM     } from '../modules/cluster_swarm'
include { BUILD_SINTAX_DB   } from '../modules/build_sintax_db'
include { TAXONOMY_SINTAX   } from '../modules/taxonomy_sintax'
include { COLLECT_RESULTS   } from '../modules/collect_results'
include { TEST_LENGTH_WF    } from './test_length'

// =============================================================================
// MAIN CECILIA WORKFLOW
// =============================================================================
workflow CECILIA {

    // -------------------------------------------------------------------------
    // INPUT CHANNEL
    // Paired-end: fromFilePairs matches *_R1_* / *_R2_* or *_1.fastq / *_2.fastq
    // Single-end: fromPath with a per-sample tuple
    // -------------------------------------------------------------------------
    if (params.paired) {
        ch_raw = Channel
            .fromFilePairs("${params.rawdata}/*_{R1,R2}*.fastq{,.gz,.bz2}", checkIfExists: true)
    } else {
        ch_raw = Channel
            .fromPath("${params.rawdata}/*.fastq{,.gz,.bz2}", checkIfExists: true)
            .map { f -> [ f.simpleName, [f] ] }
    }

    // -------------------------------------------------------------------------
    // STEP 01  —  Decompress  (per sample, parallel)
    // -------------------------------------------------------------------------
    DECOMPRESS(ch_raw)

    // -------------------------------------------------------------------------
    // STEP 02  —  FastQC on raw reads  (parallel branch, does not block DAG)
    // -------------------------------------------------------------------------
    FASTQC(DECOMPRESS.out.reads)

    // -------------------------------------------------------------------------
    // STEP 03  —  PhiX removal  (per sample)
    // -------------------------------------------------------------------------
    REMOVE_PHIX(DECOMPRESS.out.reads)

    // -------------------------------------------------------------------------
    // STEP 04  —  Paired-read assembly  (optional, per sample)
    // -------------------------------------------------------------------------
    if (params.assemble) {
        ASSEMBLE_READS(REMOVE_PHIX.out.reads)
        ch_after_assemble = ASSEMBLE_READS.out.reads
    } else {
        ch_after_assemble = REMOVE_PHIX.out.reads
    }

    // -------------------------------------------------------------------------
    // STEP 05  —  Read renaming  (optional)
    // RENAME_READS is a collect step: all samples go in together so
    // sequential numbering (sample001, sample002, ...) is assigned correctly.
    // Outputs are re-parallelised with flatten().map() for downstream steps.
    // -------------------------------------------------------------------------
    if (params.rename) {
        RENAME_READS(ch_after_assemble.map { id, f -> f }.collect())
        ch_after_rename = RENAME_READS.out.reads
            .flatten()
            .map { f -> [ f.baseName, f ] }
    } else {
        ch_after_rename = ch_after_assemble
    }

    // -------------------------------------------------------------------------
    // STEP 06  —  Primer/adapter stripping  (per sample)
    // then pool all stripped FASTQs into one pooled.fastq
    // -------------------------------------------------------------------------
    STRIP_PRIMERS(ch_after_rename)

    POOL_READS(STRIP_PRIMERS.out.fastq.map { id, f -> f }.collect())

    // Fork the pooled channel: it is reused much later by OTU table generation
    ch_pooled = POOL_READS.out.pooled_fastq

    // -------------------------------------------------------------------------
    // STEP 07  —  EE stats  (pooled)
    // Optionally strip left bases; output may be trimmed.fastq or pooled.fastq
    // -------------------------------------------------------------------------
    EE_STATS(ch_pooled)

    // Downstream source for filtering: trimmed.fastq if stripleft, else pooled.fastq
    // EE_STATS.out.trimmed is optional; using ch_pooled when stripleft=false avoids
    // an empty-channel problem downstream.
    ch_source = params.stripleft ? EE_STATS.out.trimmed : ch_pooled

    // -------------------------------------------------------------------------
    // STEP 08  —  Test read length  (optional subworkflow)
    // -------------------------------------------------------------------------
    if (params.test_length) {
        TEST_LENGTH_WF(ch_source)
    }

    // -------------------------------------------------------------------------
    // STEP 09  —  Filter by max expected errors  (per length, parallel)
    // All lengths (fixed list + user_length) deduplicated and run concurrently.
    // FILTER_MAXEE declares its tuple output optional: processes that produce
    // no reads simply emit nothing — no downstream failure.
    // -------------------------------------------------------------------------
    def all_lengths = (params.filter_lengths + [params.user_length]).unique()
    ch_lengths = Channel.fromList(all_lengths)

    FILTER_MAXEE(ch_lengths.combine(ch_source))

    // -------------------------------------------------------------------------
    // STEP 10  —  Dereplicate  (per filtered file, parallel)
    // Emits both uniques.fasta (for UPARSE/UNOISE/closed-ref) and
    // *_linear.fasta (for SWARM)
    // -------------------------------------------------------------------------
    // Only non-empty filter outputs reach DEREPLICATE (optional: true above)
    DEREPLICATE(FILTER_MAXEE.out.filtered)

    // Attach the single pooled.fastq to each (len, uniques.fasta) pair
    // so OTU-table steps have everything they need
    ch_uniques_pool = DEREPLICATE.out.uniques.combine(ch_pooled)

    // -------------------------------------------------------------------------
    // STEPS 11–15  —  Clustering  (all optional, all parallel)
    // -------------------------------------------------------------------------
    ch_cluster_fastas = Channel.empty()

    if (params.cluster_otu) {
        CLUSTER_UPARSE(ch_uniques_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_UPARSE.out.fasta)
    }

    if (params.cluster_asv) {
        CLUSTER_UNOISE(ch_uniques_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_UNOISE.out.asv)

        if (params.cluster_asv_to_otu) {
            CLUSTER_ASV_TO_OTU(
                CLUSTER_UNOISE.out.asv
                    .map { len, fasta -> tuple(len, fasta) }
                    .combine(ch_pooled)
            )
        }
    }

    if (params.closed_ref_otu) {
        CLOSED_REF_OTU(DEREPLICATE.out.uniques)
    }

    if (params.cluster_swarm) {
        ch_linear_pool = DEREPLICATE.out.linear.combine(ch_pooled)
        CLUSTER_SWARM(ch_linear_pool)
        ch_cluster_fastas = ch_cluster_fastas.mix(CLUSTER_SWARM.out.fasta)
    }

    // -------------------------------------------------------------------------
    // STEP 17  —  SINTAX taxonomy  (optional)
    // BUILD_SINTAX_DB runs once; TAXONOMY_SINTAX fans out over all cluster
    // FASTAs via .combine() — one SLURM job per FASTA, all parallel.
    // -------------------------------------------------------------------------
    if (params.sintax_taxonomy) {
        BUILD_SINTAX_DB(file(params.sintax_db))
        TAXONOMY_SINTAX(
            ch_cluster_fastas.combine(BUILD_SINTAX_DB.out.udb)
        )
    }

    // -------------------------------------------------------------------------
    // STEP 18  —  Collect and package results
    // Assemble two flat channels:
    //   ch_misc    — QC, EE stats, subset FASTA, file mapping, testLength
    //   ch_cluster — all cluster FASTAs, OTU tables, taxonomy outputs
    // Both are .collect()ed so COLLECT_RESULTS waits for the entire pipeline.
    // -------------------------------------------------------------------------

    // Misc results (always present)
    ch_misc = Channel.empty()
        .mix(FASTQC.out.html)
        .mix(POOL_READS.out.subset_fasta)
        .mix(EE_STATS.out.eestats)

    if (params.rename) {
        ch_misc = ch_misc.mix(RENAME_READS.out.mapping)
    }
    if (params.test_length) {
        ch_misc = ch_misc
            .mix(TEST_LENGTH_WF.out.results)
            .mix(TEST_LENGTH_WF.out.plots)
    }
    ch_misc = ch_misc.mix(FILTER_MAXEE.out.fastqc_html)

    // Cluster results (map tuples to plain paths)
    ch_cluster = Channel.empty()
        .mix(DEREPLICATE.out.uniques.map { len, f -> f })

    if (params.cluster_otu) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_UPARSE.out.fasta.map { len, f -> f })
            .mix(CLUSTER_UPARSE.out.table)
    }
    if (params.cluster_asv) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_UNOISE.out.asv.map { len, f -> f })
            .mix(CLUSTER_UNOISE.out.table)
    }
    if (params.cluster_asv_to_otu) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_ASV_TO_OTU.out.fasta.map { len, f -> f })
            .mix(CLUSTER_ASV_TO_OTU.out.table)
    }
    if (params.closed_ref_otu) {
        ch_cluster = ch_cluster
            .mix(CLOSED_REF_OTU.out.table.map { len, f -> f })
    }
    if (params.cluster_swarm) {
        ch_cluster = ch_cluster
            .mix(CLUSTER_SWARM.out.fasta.map { len, f -> f })
            .mix(CLUSTER_SWARM.out.table)
    }
    if (params.sintax_taxonomy) {
        ch_cluster = ch_cluster
            .mix(TAXONOMY_SINTAX.out.sintax.map { len, f -> f })
    }

    COLLECT_RESULTS(
        ch_misc.collect(),
        ch_cluster.collect()
    )
}
