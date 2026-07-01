/*
 * ANNOTATE
 *
 * Merges all scored FlashFry chunks back together, maps each guide back
 * to its source region (via region_metadata from PREP_REGIONS), and
 * computes each guide's distance to the target region's reference point
 * (e.g. region center / feature of interest). Guides are flagged
 * pass/fail against `params.blast_radius_bp` (the ~1kb CRISPRi silencing
 * window) — a guide's silencing effect is only guaranteed to cover the
 * intended target if it falls within that window.
 */

process ANNOTATE {
    label 'utils_process'
    tag "merge"

    container 'TODO/crispri-guide-designer-utils:tag'

    input:
    path scored_chunks   // collected list of all scored.output files
    path region_metadata
    val  blast_radius_bp

    output:
    path "guides.annotated.tsv", emit: annotated_guides

    script:
    """
    annotate.py \\
        --scored_chunks ${scored_chunks} \\
        --region_metadata ${region_metadata} \\
        --blast_radius_bp ${blast_radius_bp} \\
        --out guides.annotated.tsv
    """
}
