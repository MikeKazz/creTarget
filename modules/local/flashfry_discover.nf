/*
 * FLASHFRY_DISCOVER
 *
 * Runs FlashFry's `discover` mode against the cached hg38 index to find
 * all candidate NGG-PAM guides (+ their off-target hits) within each
 * region's flanking FASTA. PAM-finding is delegated entirely to FlashFry
 * here (per design decision) rather than the R helper in the original
 * script.
 *
 * Nextflow scatters this automatically across however many region-FASTA
 * chunks are emitted upstream — no fixed chunk count needed (replaces the
 * manual SLURM --array=0-28 chunking approach).
 */

process FLASHFRY_DISCOVER {
    label 'flashfry_discover'
    tag "${region_fasta.simpleName}"

    container 'TODO/flashfry-image:1.15'   // custom image with FlashFry-assembly-1.15.jar baked in; ignored under -profile bare

    input:
    path region_fasta
    path flashfry_db
    val  max_mismatch
    val  flashfry_jar   // path string; container image's baked-in JAR used if this equals the default container path

    output:
    path "discover.output", emit: discover_output

    script:
    """
    java -Xmx14g -jar ${flashfry_jar} \\
        discover \\
        --database ${flashfry_db} \\
        --fasta ${region_fasta} \\
        --output discover.output \\
        --maxMismatch ${max_mismatch}
    """
}
