/*
 * PREP_REGIONS
 *
 * Input: user's BED file of target regions + reference FASTA.
 * - Optionally collapses targets within `params.blast_radius_bp` of each
 *   other (CRISPRi silencing radius), keeping one representative per
 *   cluster and recording which targets were treated as redundant.
 *   (Ports the iterative dedup logic from the original hFOB SNP script,
 *   generalized to arbitrary BED regions/any "closeness" metric provided
 *   in the BED, rather than ATAC-seq distance specifically.)
 * - Extracts a per-region FASTA with `params.flank_bp` of flanking sequence
 *   on each side (flanking context is required by FlashFry's on-target
 *   scoring metrics), ready to hand to FLASHFRY_DISCOVER.
 *
 * Output FASTA headers encode region_id + original coordinates so
 * downstream steps can map guides back to their source region.
 */

process PREP_REGIONS {
    label 'utils_process'
    tag "${bed.simpleName}"

    container 'TODO/crispri-guide-designer-utils:tag'

    input:
    tuple val(input_type), path(bed)   // input_type: 'snp' or 'enhancer'
    path fasta
    path fasta_index
    val  dedup_blast_radius
    val  blast_radius_bp
    val  flank_bp
    val  snp_window_bp

    output:
    path "regions.fasta",        emit: region_fasta
    path "regions.metadata.tsv", emit: region_metadata

    script:
    def dedup_flag = dedup_blast_radius ? "--dedup_blast_radius" : ""
    """
    prep_regions.py \\
        --bed ${bed} \\
        --input_type ${input_type} \\
        --fasta ${fasta} \\
        --flank_bp ${flank_bp} \\
        --snp_window_bp ${snp_window_bp} \\
        --link_distance_bp ${blast_radius_bp} \\
        ${dedup_flag} \\
        --outdir .
    """
}
