#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * crispri-guide-designer
 *
 * Usage (once reference bundle + Docker images are published):
 *   nextflow run main.nf --snp_bed snps.bed --enhancer_bed enhancers.bed \
 *       [-profile standard|slurm] [-profile docker|singularity]
 *
 * Usage (testing now, against an existing local FlashFry install):
 *   nextflow run main.nf --snp_bed snps.bed --enhancer_bed enhancers.bed \
 *       --local_fasta /path/to/hg38.fa \
 *       --local_flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
 *       --local_flashfry_db /path/to/hg38_standard_cas9ngg_db \
 *       -profile bare
 *
 * At least one of --snp_bed / --enhancer_bed is required; both may be
 * supplied together. If any --local_* param is set, SETUP_REFERENCE is
 * skipped entirely and those paths are used directly.
 */

include { SETUP_REFERENCE   } from './modules/local/setup_reference.nf'
include { PREP_REGIONS as PREP_SNP_REGIONS      } from './modules/local/prep_regions.nf'
include { PREP_REGIONS as PREP_ENHANCER_REGIONS } from './modules/local/prep_regions.nf'
include { FLASHFRY_DISCOVER } from './modules/local/flashfry_discover.nf'
include { FLASHFRY_SCORE    } from './modules/local/flashfry_score.nf'
include { ANNOTATE          } from './modules/local/annotate.nf'
include { RANK_AND_REPORT   } from './modules/local/rank_report.nf'

def helpMessage() {
    log.info """
    crispri-guide-designer
    -----------------------
    Usage:
      nextflow run main.nf [--snp_bed <snps.bed>] [--enhancer_bed <enhancers.bed>] [-profile standard|slurm|bare]

    At least one of --snp_bed / --enhancer_bed is required.

    Optional:
      --dedup_blast_radius   Flag targets within blast_radius_bp of each other (does not merge) [default: false]
      --blast_radius_bp      CRISPRi silencing window, bp [default: 1000]
      --flank_bp              Enhancer-type: padding added around each BED interval [default: 50]
      --snp_window_bp          SNP-type: search window on each side of the point [default: 150]
      --cache_dir             Where the reference genome + FlashFry index are cached
                               [default: ~/.crispri-guide-designer]
      --outdir                 Output directory [default: results]

    Local-install testing (bypasses SETUP_REFERENCE, use with -profile bare):
      --local_fasta            Path to an existing hg38 FASTA
      --local_flashfry_jar     Path to an existing FlashFry-assembly-*.jar
      --local_flashfry_db      Path to an existing FlashFry index directory
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.snp_bed && !params.enhancer_bed) {
        error "At least one of --snp_bed or --enhancer_bed is required. Run with --help for usage."
    }

    def use_local = params.local_fasta || params.local_flashfry_jar || params.local_flashfry_db
    if (use_local && !(params.local_fasta && params.local_flashfry_jar && params.local_flashfry_db)) {
        error "If any --local_* param is set, all three (--local_fasta, --local_flashfry_jar, --local_flashfry_db) must be set."
    }

    // ---- 1. Reference + FlashFry index: either local override or auto-cached ----
    def fasta_ch
    def fasta_index_ch
    def flashfry_db_ch
    def flashfry_jar_val

    if (use_local) {
        fasta_ch         = Channel.fromPath(params.local_fasta, checkIfExists: true)
        fasta_index_ch    = Channel.fromPath("${params.local_fasta}.fai", checkIfExists: true)
        flashfry_db_ch     = Channel.fromPath(params.local_flashfry_db, checkIfExists: true)
        flashfry_jar_val    = params.local_flashfry_jar
    } else {
        SETUP_REFERENCE(params.reference_version)
        fasta_ch         = SETUP_REFERENCE.out.fasta
        fasta_index_ch    = SETUP_REFERENCE.out.fasta_index
        flashfry_db_ch     = SETUP_REFERENCE.out.flashfry_db
        flashfry_jar_val    = '/opt/flashfry/FlashFry-assembly-1.15.jar'   // baked into the container image
    }

    // ---- 2. Parse BED(s), optional link-flagging, extract per-region flanking FASTA ----
    region_fasta_ch    = Channel.empty()
    region_metadata_ch = Channel.empty()

    if (params.snp_bed) {
        PREP_SNP_REGIONS(
            Channel.of(['snp', file(params.snp_bed, checkIfExists: true)]),
            fasta_ch,
            fasta_index_ch,
            params.dedup_blast_radius,
            params.blast_radius_bp,
            params.flank_bp,
            params.snp_window_bp
        )
        region_fasta_ch    = region_fasta_ch.mix(PREP_SNP_REGIONS.out.region_fasta)
        region_metadata_ch = region_metadata_ch.mix(PREP_SNP_REGIONS.out.region_metadata)
    }

    if (params.enhancer_bed) {
        PREP_ENHANCER_REGIONS(
            Channel.of(['enhancer', file(params.enhancer_bed, checkIfExists: true)]),
            fasta_ch,
            fasta_index_ch,
            params.dedup_blast_radius,
            params.blast_radius_bp,
            params.flank_bp,
            params.snp_window_bp
        )
        region_fasta_ch    = region_fasta_ch.mix(PREP_ENHANCER_REGIONS.out.region_fasta)
        region_metadata_ch = region_metadata_ch.mix(PREP_ENHANCER_REGIONS.out.region_metadata)
    }

    // ---- 3. FlashFry discover (candidate guides + raw off-target hits) ----
    // TODO: scatter across multiple FASTA chunks for large inputs rather than
    // running discover once per input-type channel element.
    FLASHFRY_DISCOVER(
        region_fasta_ch,
        flashfry_db_ch,
        params.flashfry_max_mismatch,
        flashfry_jar_val
    )

    // ---- 4. FlashFry score (on-target + off-target metrics) ----
    FLASHFRY_SCORE(
        FLASHFRY_DISCOVER.out.discover_output,
        flashfry_db_ch,
        params.flashfry_max_mismatch,
        params.flashfry_scoring_metrics,
        flashfry_jar_val
    )

    // ---- 5. Merge + annotate against source regions ----
    // Genome coordinates are recovered by parsing the `contig` field
    // (region_id__chrom:start-end) written into the FASTA header by
    // PREP_REGIONS -- confirmed working against real FlashFry output.
    ANNOTATE(
        FLASHFRY_SCORE.out.scored_output.collect(),
        region_metadata_ch.collectFile(name: 'regions.metadata.merged.tsv', keepHeader: true),
        params.blast_radius_bp
    )

    // ---- 6. Rank + final report ----
    RANK_AND_REPORT(
        ANNOTATE.out.annotated_guides,
        params.flashfry_max_mismatch
    )

    RANK_AND_REPORT.out.ranked_guides
        .subscribe { f -> log.info "Final ranked guides written to: ${f}" }
}
