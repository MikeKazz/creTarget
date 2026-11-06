/*
 * Writes a provenance record next to the results.
 *
 * Every output folder should be able to explain itself a year later: which
 * version of the pipeline produced it, with which parameters, using which
 * tool versions. Reviewers ask; so does future-you.
 *
 * This is a `script:` process rather than `exec:` because exec blocks run
 * inside the Nextflow JVM and do not stage input files, so the collected
 * versions_*.yml would not exist on disk.
 */
process RUN_MANIFEST {
    tag "manifest"
    label 'process_single'

    input:
    // Every process emits a file literally named versions.yml, so they must be
    // staged under distinct names or Nextflow rejects the collision.
    path versions, stageAs: 'versions_*.yml'

    output:
    path "run_manifest.json"   , emit: manifest
    path "software_versions.yml", emit: software_versions

    script:
    def manifest_json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        pipeline: [
            name     : workflow.manifest.name,
            version  : workflow.manifest.version,
            revision : workflow.revision ?: 'local',
            commit_id: workflow.commitId ?: 'n/a',
            home_page: workflow.manifest.homePage,
        ],
        run: [
            run_name        : workflow.runName,
            session_id      : workflow.sessionId.toString(),
            start_time      : workflow.start.toString(),
            command_line    : workflow.commandLine,
            profile         : workflow.profile,
            container_engine: workflow.containerEngine ?: 'none',
            nextflow_version: workflow.nextflow.version.toString(),
            work_dir        : workflow.workDir.toString(),
            launch_dir      : workflow.launchDir.toString(),
            user            : workflow.userName,
        ],
        params: params.findAll { k, v -> v != null && k != 'genomes' },
    ]))
    """
    cat <<'END_MANIFEST' > run_manifest.json
${manifest_json}
END_MANIFEST

    # Collapse the per-process versions files into one document. Deduplication
    # is per BLOCK, keyed on the process-name header -- one versions.yml comes
    # back per chunk, so the same block repeats N times. Deduplicating by line
    # instead would drop a later process's "cretarget: 0.2.0" just because an
    # earlier block already had that exact line, leaving empty sections.
    cat versions_*.yml | awk '
        /^"/ { header = \$0; skip = (header in seen); seen[header] = 1 }
        !skip { print }
    ' > software_versions.yml
    """
}
