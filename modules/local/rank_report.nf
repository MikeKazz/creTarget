/*
 * RANK_AND_REPORT
 *
 * Implements the project owner's real ranking scheme (see bin/rank_report.py
 * for full detail): hard filters (seed_score, dangerous_GC, dangerous_polyT,
 * multi-site genomic match, FlashFry overflow) mark guides as
 * `exclude (<reason>)` in the `rank` column; guides passing all filters get
 * a locally-ranked (per region_id) weighted composite score.
 *
 * NOTE: the intended scheme also weights "launch-dCas9" on-target efficacy
 * scores (42% of total weight, for both enhancers and promoters) -- that
 * tool isn't wired into this pipeline yet, so it's omitted for now and the
 * remaining weights are renormalized. Ranks will change once it's added.
 */

process RANK_AND_REPORT {
    label 'utils_process'
    tag "final"

    container 'TODO/crispri-guide-designer-utils:tag'

    input:
    path annotated_guides
    val  max_mismatch

    output:
    path "guides.ranked.tsv", emit: ranked_guides

    script:
    """
    rank_report.py \\
        --annotated ${annotated_guides} \\
        --max_mismatch ${max_mismatch} \\
        --out guides.ranked.tsv
    """
}
