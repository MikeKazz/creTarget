/*
 * Promoter-proximal seed screening.
 *
 * Runs ONCE over the collated table, not per chunk. The scan's cost is one
 * pass over the promoter reference and barely moves with guide count, so
 * per-chunk invocation paid that fixed pass on every chunk -- at chunk_size
 * 1000 a genome-wide run would have spent weeks re-reading the same ~212 Mb.
 *
 * `--seed_finder_max_guides` splits the table across several tasks. The limit
 * it exists for is memory, not time: the seed index is the only thing held in
 * RAM, at roughly 1 GB per million guides screened.
 */
process SEED_FINDER {
    tag "${prefix}"

    input:
    path ranked
    path promoter_fa
    path tss_table

    output:
    path "${prefix}.final.tsv", emit: final_table
    path "${prefix}_raw_hits.csv", emit: raw_hits
    path "versions.yml", emit: versions

    script:
    prefix = ranked.name.replaceAll(/\.(tsv|txt|csv)$/, '')
    def all_guides = params.seed_finder_all_guides ? '--all-guides' : ''
    """
    cretarget seed-finder \\
        --ranked ${ranked} \\
        --promoter-fa ${promoter_fa} \\
        --tss-table ${tss_table} \\
        --min-seed-length ${params.seed_finder_min_length} \\
        --tss-window ${params.seed_finder_tss_window} \\
        ${all_guides} \\
        --outdir seed_finder_work \\
        --raw-hits ${prefix}_raw_hits.csv \\
        --out ${prefix}.final.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
        seed_reference: "\$(test -f ${tss_table} && wc -l < ${tss_table} | tr -d ' ' || echo '?') TSS rows"
    END_VERSIONS
    """
}
