/*
 * seed-finder-lite: PAM-proximal seed / TSS-proximity screening.
 *
 * This runs ONCE over the whole ranked table, not once per chunk. The tool's
 * batch mode scans the promoterome a single time and matches all guides
 * against it by dictionary lookup, so its runtime is effectively independent
 * of guide count (~25 s for 20 guides, ~25 s for 2,000). Per-chunk invocation
 * paid that fixed scan on every chunk, which at chunk_size 1000 meant weeks of
 * re-reading the same 212 MB for a genome-wide run.
 *
 * `--seed_finder_max_guides` splits the table into several tasks for cases
 * where one task's memory or wall time is a problem; each shard pays its own
 * scan, so keep the shard count low.
 */
process SEED_FINDER {
    tag "${prefix}"
    label 'process_low'

    input:
    path ranked
    path seed_finder_dir

    output:
    path "${prefix}.final.tsv", emit: final_table
    path "${prefix}_raw_hits.csv", emit: raw_hits
    path "versions.yml", emit: versions

    script:
    prefix = ranked.name.replaceAll(/\.(tsv|txt|csv)$/, '')
    def python_arg = params.seed_finder_python ? "--python-bin ${params.seed_finder_python}" : ''
    /*
     * The self-test asks the tool to recover two real Cas9 sites taken from
     * its own promoter reference, one per orientation. It exists because a
     * build whose plus-strand scan indexes guides under the reverse
     * complement of their seed reports only minus-strand hits -- it exits 0
     * and returns plausible numbers while missing about half of all true
     * hits. One extra promoterome scan is cheap insurance against silently
     * publishing those results.
     */
    // Kept on one line on purpose: a line-continuation backslash inside this
    // interpolated string needs double-escaping in Groovy, and getting that
    // wrong emits a literal "\\" that bash reads as an escaped backslash
    // rather than a continuation.
    def preflight = params.seed_finder_skip_check ? '' :
        "cretarget seed-finder-check --seed-finder-dir ${seed_finder_dir} ${python_arg} --outdir seed_finder_check"
    """
    ${preflight}
    cretarget seed-finder \\
        --ranked ${ranked} \\
        --seed-finder-dir ${seed_finder_dir} \\
        --min-seed-length ${params.seed_finder_min_length} \\
        --tss-window ${params.seed_finder_tss_window} \\
        ${python_arg} \\
        --outdir seed_finder_work \\
        --raw-hits ${prefix}_raw_hits.csv \\
        --out ${prefix}.final.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
        seed-finder-lite: "\$(cat ${seed_finder_dir}/VERSION 2>/dev/null || echo 'unknown (no VERSION file)')"
    END_VERSIONS
    """
}
