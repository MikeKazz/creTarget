#!/usr/bin/env nextflow
/*
 * creTarget
 * Design, score, and rank CRISPRi guide RNAs for noncoding target regions.
 *
 * https://github.com/MikeKazz/creTarget
 */

nextflow.enable.dsl = 2

include { CRETARGET } from './workflows/cretarget'

def helpMessage() {
    log.info """
    creTarget ${workflow.manifest.version}

    Usage:
      nextflow run MikeKazz/creTarget -profile <local|slurm|pbs|sge|lsf>[,conda] \\
          --enhancer_bed enhancers.txt \\
          --fasta hg38.fa \\
          --flashfry_jar FlashFry-assembly-1.15.jar \\
          --flashfry_db hg38_standard_cas9ngg_db \\
          --outdir results

    Input (at least one required):
      --snp_bed                 SNP targets: 'chr:pos' per line (1-based), or BED
      --enhancer_bed            Enhancer targets: 'chr:start-end' per line, or BED

    Reference (required):
      --fasta                   hg38 FASTA (a .fai alongside it is expected)
      --flashfry_jar            FlashFry assembly JAR
      --flashfry_db             FlashFry off-target index (its .header must sit alongside)

    Region prep:
      --flank_bp                Padding around each enhancer interval  [${params.flank_bp}]
      --snp_window_bp           Search window each side of a SNP       [${params.snp_window_bp}]
      --blast_radius_bp         CRISPRi silencing window               [${params.blast_radius_bp}]
      --dedup_blast_radius      Flag targets within one blast radius of each other

    Scale:
      --chunk_size              Regions per parallel task              [${params.chunk_size}]

    Optional annotation:
      --atac_narrowpeak         narrowPeak file; adds distance_to_atac_summit
      --run_launch_dcas9        Add launch-dCas9 CNN on-target scoring
        --launch_script           Path to launch-dCas9.py
        --launch_python           Python from a conda env with torch
        --ckpt_dir                Checkpoint dir (see `cretarget fetch-models`)
      --run_seed_finder         Add seed-finder-lite off-target/TSS columns
        --seed_finder_dir         seed-finder-lite directory
        --seed_finder_min_length  Minimum exact-match seed length (no default)
        --seed_finder_tss_window  Max bp from a TSS to count as a hit (no default)

    Output:
      --outdir                  Output directory                       [${params.outdir}]

    Nextflow flags worth knowing: -resume (reuse completed tasks),
    -profile, -work-dir, -with-report.
    """.stripIndent()
}

def validateParams() {
    def errors = []

    if (!params.snp_bed && !params.enhancer_bed) {
        errors << "At least one of --snp_bed / --enhancer_bed is required."
    }
    if (!params.fasta)        errors << "--fasta is required."
    if (!params.flashfry_jar) errors << "--flashfry_jar is required."
    if (!params.flashfry_db)  errors << "--flashfry_db is required."

    if (params.fasta && !file("${params.fasta}.fai").exists()) {
        errors << "No FASTA index at ${params.fasta}.fai -- run: samtools faidx ${params.fasta}"
    }
    if (params.flashfry_db && !file("${params.flashfry_db}.header").exists()) {
        errors << "FlashFry index header not found at ${params.flashfry_db}.header " +
                  "(it must sit alongside the database file)."
    }

    if (params.run_launch_dcas9) {
        if (!params.launch_script) errors << "--run_launch_dcas9 requires --launch_script."
        if (!params.ckpt_dir) {
            errors << "--run_launch_dcas9 requires --ckpt_dir. Fetch the checkpoints " +
                      "first with: cretarget fetch-models"
        }
    }

    if (params.run_seed_finder) {
        if (!params.seed_finder_dir) {
            errors << "--run_seed_finder requires --seed_finder_dir."
        }
        // Deliberately no defaults: these change what the screen means, so
        // the user has to choose them consciously.
        if (params.seed_finder_min_length == null) {
            errors << "--run_seed_finder requires --seed_finder_min_length (no default; " +
                      "see seed-finder-lite's documentation)."
        }
        if (params.seed_finder_tss_window == null) {
            errors << "--run_seed_finder requires --seed_finder_tss_window (no default)."
        }
    }

    if (params.chunk_size != null && (params.chunk_size as int) < 1) {
        errors << "--chunk_size must be >= 1."
    }

    /*
     * A container engine with no image is the dangerous case: Nextflow would
     * run every task uncontainerised, using whatever tools happen to be on
     * the host, while the run reports a container engine and looks
     * reproducible. Refuse rather than produce quietly untrustworthy results.
     */
    if (workflow.containerEngine && !params.container_image) {
        errors << "Profile '${workflow.profile}' enables the " +
                  "${workflow.containerEngine} container engine but no image was " +
                  "given. Pass --container_image <image>, or drop the container " +
                  "profile and use -profile conda (or a host install) instead. " +
                  "Running with an engine enabled and no image would silently " +
                  "execute tasks against host tools."
    }
    if (params.container_image && !workflow.containerEngine) {
        log.warn "--container_image was given (${params.container_image}) but no " +
                 "container engine is enabled; it will be ignored. Add " +
                 "-profile singularity (or docker/apptainer)."
    }

    if (errors) {
        log.error "Parameter validation failed:\n" + errors.collect { "  - ${it}" }.join('\n')
        log.error "\nRun with --help for usage."
        System.exit(1)
    }
}

workflow {

    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    if (params.version) {
        log.info "creTarget ${workflow.manifest.version}"
        System.exit(0)
    }

    validateParams()

    log.info """
    ================================================================
     creTarget ${workflow.manifest.version}
    ================================================================
     SNP targets      : ${params.snp_bed ?: '-'}
     Enhancer targets : ${params.enhancer_bed ?: '-'}
     Reference        : ${params.fasta}
     Chunk size       : ${params.chunk_size} regions/task
     launch-dCas9     : ${params.run_launch_dcas9 ? 'enabled' : 'disabled'}
     seed-finder-lite : ${params.run_seed_finder ? 'enabled' : 'disabled'}
     Executor         : ${workflow.profile}
     Output           : ${params.outdir}
    ================================================================
    """.stripIndent()

    /*
     * Optional file inputs are passed as named placeholders rather than
     * empty channels: a process input path is mandatory in DSL2, and a
     * placeholder keeps the process signature uniform while letting the
     * script decide whether the argument was really supplied.
     */
    ch_snp_bed = params.snp_bed
        ? Channel.fromPath(params.snp_bed, checkIfExists: true).first()
        : Channel.fromPath("${projectDir}/assets/NO_SNP_FILE").first()

    ch_enhancer_bed = params.enhancer_bed
        ? Channel.fromPath(params.enhancer_bed, checkIfExists: true).first()
        : Channel.fromPath("${projectDir}/assets/NO_ENH_FILE").first()

    ch_atac = params.atac_narrowpeak
        ? Channel.fromPath(params.atac_narrowpeak, checkIfExists: true).first()
        : Channel.fromPath("${projectDir}/assets/NO_ATAC_FILE").first()

    ch_fasta = Channel.of(
        [file(params.fasta, checkIfExists: true), file("${params.fasta}.fai", checkIfExists: true)]
    ).first()

    ch_flashfry_jar = Channel.fromPath(params.flashfry_jar, checkIfExists: true).first()

    ch_flashfry_db = Channel.of([
        file(params.flashfry_db, checkIfExists: true),
        file("${params.flashfry_db}.header", checkIfExists: true)
    ]).first()

    CRETARGET(
        ch_snp_bed,
        ch_enhancer_bed,
        ch_fasta,
        ch_flashfry_jar,
        ch_flashfry_db,
        ch_atac
    )
}
