/*
 * creTarget main workflow.
 *
 * Shape of the thing:
 *
 *   PREP_REGIONS (once, global -- linked-region flagging needs every target)
 *        |
 *   splitFasta  ->  N chunks
 *        |
 *   FLASHFRY_DISCOVER_SCORE  \
 *   ANNOTATE                  |  all per-chunk, fully parallel,
 *   LAUNCH_DCAS9 (optional)   |  no cross-chunk synchronization
 *   RANK                      |
 *   SEED_FINDER (optional)   /
 *        |
 *   collectFile -> one output table
 *
 * Everything downstream of splitting is scoped to a single region_id, and
 * every guide for a region lands in the same chunk, so ranking is correct
 * per-chunk and no global merge is needed. That is what makes memory flat
 * with respect to target count: 3,000 regions and 2.35M regions use the
 * same per-task memory, only the task count differs.
 */

include { PREP_REGIONS            } from '../modules/local/prep_regions'
include { FLASHFRY_DISCOVER_SCORE } from '../modules/local/flashfry'
include { ANNOTATE                } from '../modules/local/annotate'
include { LAUNCH_DCAS9            } from '../modules/local/launch_dcas9'
include { RANK                    } from '../modules/local/rank'
include { SEED_FINDER             } from '../modules/local/seed_finder'
include { RUN_MANIFEST            } from '../modules/local/run_manifest'

workflow CRETARGET {

    take:
    ch_snp_bed
    ch_enhancer_bed
    ch_fasta          // tuple(fasta, fai)
    ch_flashfry_jar
    ch_flashfry_db    // tuple(db, db.header)
    ch_atac

    main:
    ch_versions = Channel.empty()

    //
    // Prepare regions once, globally.
    //
    PREP_REGIONS(ch_snp_bed, ch_enhancer_bed, ch_fasta)
    ch_versions = ch_versions.mix(PREP_REGIONS.out.versions)

    //
    // Split into chunks. splitFasta streams the file rather than holding
    // every record in memory, which the old Python splitter did not.
    //
    ch_chunks = PREP_REGIONS.out.fasta
        .splitFasta(by: params.chunk_size, file: true)

    //
    // FlashFry discover + score, one task per chunk.
    //
    FLASHFRY_DISCOVER_SCORE(ch_chunks, ch_flashfry_jar, ch_flashfry_db)
    ch_versions = ch_versions.mix(FLASHFRY_DISCOVER_SCORE.out.versions.first())

    //
    // Annotate per chunk. The metadata table is global; it is consumed as a
    // value channel so every task sees it, and the step streams it rather
    // than loading it whole.
    //
    ANNOTATE(
        FLASHFRY_DISCOVER_SCORE.out.scored,
        PREP_REGIONS.out.metadata,
        ch_atac
    )
    ch_versions = ch_versions.mix(ANNOTATE.out.versions.first())

    ch_for_rank = ANNOTATE.out.annotated

    //
    // Optional launch-dCas9 CNN on-target scoring.
    //
    if (params.run_launch_dcas9) {
        LAUNCH_DCAS9(
            ANNOTATE.out.annotated,
            Channel.fromPath(params.launch_script, checkIfExists: true).first(),
            Channel.fromPath(params.ckpt_dir, type: 'dir', checkIfExists: true).first()
        )
        ch_for_rank = LAUNCH_DCAS9.out.annotated
        ch_versions = ch_versions.mix(LAUNCH_DCAS9.out.versions.first())
    }

    //
    // Hard filters + composite score + per-region rank.
    //
    RANK(ch_for_rank)
    ch_versions = ch_versions.mix(RANK.out.versions.first())

    ch_final = RANK.out.ranked

    //
    // Optional seed-finder-lite off-target / TSS screening.
    //
    if (params.run_seed_finder) {
        SEED_FINDER(
            RANK.out.ranked,
            Channel.fromPath(params.seed_finder_dir, type: 'dir', checkIfExists: true).first()
        )
        ch_final = SEED_FINDER.out.final_table
        ch_versions = ch_versions.mix(SEED_FINDER.out.versions.first())
    }

    //
    // One table out. collectFile streams the per-chunk files together and
    // keeps a single header, so this stays cheap even with thousands of
    // chunks -- nothing is ever all in memory at once.
    //
    def output_name = params.run_seed_finder ? 'guides.final.tsv' : 'guides.ranked.tsv'

    ch_collected = ch_final
        .map { chunk_id, table -> table }
        .collectFile(
            name: output_name,
            storeDir: params.outdir,
            keepHeader: true,
            skip: 1,
            /*
             * Sort by the chunk's numeric index, NOT by `sort: true`.
             *
             * `sort: true` orders by the Path's natural ordering, which is the
             * absolute path -- and that contains the task's random work-directory
             * hash. Two identical runs therefore concatenated the chunks in
             * different orders and produced byte-different output files. No rows
             * were lost and every row is self-labelled, but "same input, same
             * output" is a claim this pipeline has to be able to make.
             *
             * Per-chunk files are named "<prefix>.<n>.ranked.tsv". Zero-padding
             * the index makes the string order match the numeric order, so chunk
             * 2 precedes chunk 10 rather than following chunk 1. The fallback
             * also returns a String, so the comparator never mixes types.
             */
            sort: { f ->
                def m = (f.name =~ /\.(\d+)\./)
                m.find() ? String.format('%09d', m.group(1) as Integer) : f.name
            }
        )

    //
    // Provenance record.
    //
    RUN_MANIFEST(ch_versions.collect())

    emit:
    guides   = ch_collected
    manifest = RUN_MANIFEST.out.manifest
    versions = ch_versions
}
