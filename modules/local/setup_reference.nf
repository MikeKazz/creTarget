/*
 * SETUP_REFERENCE
 *
 * Ensures a hg38 reference FASTA (main chromosomes only, no alt scaffolds,
 * no chrM) and a matching FlashFry index are present in the local cache
 * directory. Downloads them from the hosted release bundle on first run;
 * on subsequent runs, detects the cache and emits paths immediately without
 * re-downloading.
 *
 * NOTE: this process intentionally does NOT run inside the FlashFry
 * container — it just needs curl/wget + checksum tools, so it uses a
 * minimal utility image.
 */

process SETUP_REFERENCE {
    label 'setup_reference'
    tag "${params.reference_version}"

    container 'quay.io/biocontainers/curl:...'   // TODO: pin a minimal curl/coreutils image

    // Cache lives outside Nextflow's work dir on purpose (persists across
    // separate pipeline runs/projects, not just within one run's -resume).
    storeDir params.cache_dir

    input:
    val reference_version

    output:
    path "${reference_version}/genome.fa",          emit: fasta
    path "${reference_version}/genome.fa.fai",      emit: fasta_index
    path "${reference_version}/flashfry_db",        emit: flashfry_db
    path "${reference_version}/flashfry_db.paramfile", emit: flashfry_db_paramfile, optional: true

    script:
    """
    # TODO:
    # 1. Check if genome.fa / flashfry_db already exist under storeDir (Nextflow
    #    storeDir + emit will short-circuit reruns automatically once cached).
    # 2. If absent, download from \${params.reference_url_base}/\${reference_version}/
    #    (FASTA + prebuilt FlashFry index bundle, hosted on Zenodo/GitHub Releases).
    # 3. Verify checksums.
    # 4. Unpack into ${reference_version}/.
    echo "TODO: implement download + cache logic"
    """
}
