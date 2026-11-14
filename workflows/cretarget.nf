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
 *   RANK                      /
 *        |
 *   collectFile -> guides.ranked.tsv
 *        |
 *   SEED_FINDER (optional, ONCE over the whole table)
 *        |
 *   collectFile -> guides.final.tsv
 *
 * SEED_FINDER is deliberately NOT per-chunk. seed-finder-lite's batch mode
 * scans the promoterome once and then matches every guide by dictionary
 * lookup, so its cost barely moves with guide count (~25 s for 20 guides,
 * ~25 s for 2,000). Running it per chunk paid that scan on every chunk; at
 * chunk_size 1000 a genome-wide run would have spent weeks re-reading the
 * same 212 MB. It therefore sits after collation, where one invocation covers
 * every guide. The step streams the collated table rather than loading it, so
 * this does not trade runtime for memory.
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
    /*
     * `as int` is load-bearing. Nextflow 26 stopped coercing command-line
     * params to numbers, so `--chunk_size 1000` arrives as a String while the
     * same value from a config file or a profile arrives as an Integer.
     * splitFasta rejects the String with "Value don't match: class
     * java.lang.Integer", which reads as if Integer were the problem when it
     * is in fact what is expected. Nextflow 25 and earlier coerced it, so this
     * only breaks for people on a current Nextflow.
     */
    ch_chunks = PREP_REGIONS.out.fasta
        .splitFasta(by: params.chunk_size as int, file: true)

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

    //
    // Collate the per-chunk ranked tables. collectFile streams them together
    // and keeps a single header, so this stays cheap even with thousands of
    // chunks -- nothing is ever all in memory at once.
    //
    ch_ranked = RANK.out.ranked
        .map { chunk_id, table -> table }
        .collectFile(
            name: 'guides.ranked.tsv',
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
    // Optional seed-finder-lite off-target / TSS screening, over the collated
    // table. Sharding is opt-in: each shard repeats the promoterome scan, so
    // it is only worth it when one task cannot hold the hit summary or finish
    // inside its walltime.
    //
    ch_collected = ch_ranked

    if (params.run_seed_finder) {
        /*
         * The promoter reference is built once with `cretarget
         * build-seed-reference` and lives in the cache. Resolving it here
         * rather than inside the process means a missing reference is a
         * channel error at launch, not a task failure after the whole ranking
         * stage has already run.
         */
        def cache_root = System.getenv('CRETARGET_CACHE_DIR')
            ?: "${System.getProperty('user.home')}/.cretarget"
        ch_seed_promoter_fa = params.seed_finder_promoter_fa
            ?: "${cache_root}/seed_reference/promoter_windows.fa.gz"
        ch_seed_tss_table = params.seed_finder_tss_table
            ?: "${cache_root}/seed_reference/tss_table.tsv"

        ch_seed_input = params.seed_finder_max_guides
            ? ch_ranked.splitText(
                  by: params.seed_finder_max_guides as int,
                  keepHeader: true,
                  file: true
              )
            : ch_ranked

        SEED_FINDER(
            ch_seed_input,
            Channel.fromPath(ch_seed_promoter_fa, checkIfExists: true).first(),
            Channel.fromPath(ch_seed_tss_table, checkIfExists: true).first()
        )
        ch_versions = ch_versions.mix(SEED_FINDER.out.versions.first())

        /*
         * Collated even when there is a single shard: one code path, and the
         * published name stays the same either way. The sort mirrors the
         * ranked collation above -- numeric by shard index, never by path,
         * because a work-directory hash in the path would make output order
         * vary between otherwise identical runs.
         */
        ch_collected = SEED_FINDER.out.final_table
            .collectFile(
                name: 'guides.final.tsv',
                storeDir: params.outdir,
                keepHeader: true,
                skip: 1,
                sort: { f ->
                    def m = (f.name =~ /\.(\d+)\./)
                    m.find() ? String.format('%09d', m.group(1) as Integer) : f.name
                }
            )
    }

    //
    // Provenance record.
    //
    RUN_MANIFEST(ch_versions.collect())

    emit:
    guides   = ch_collected
    manifest = RUN_MANIFEST.out.manifest
    versions = ch_versions
}
