process SEED_FINDER {
    tag "${chunk_id}"

    input:
    tuple val(chunk_id), path(ranked)
    path seed_finder_dir

    output:
    tuple val(chunk_id), path("${chunk_id}.final.tsv"), emit: final_table
    path "${chunk_id}_raw_hits.csv", emit: raw_hits
    path "versions.yml", emit: versions

    script:
    def python_arg = params.seed_finder_python ? "--python-bin ${params.seed_finder_python}" : ''
    """
    cretarget seed-finder \\
        --ranked ${ranked} \\
        --seed-finder-dir ${seed_finder_dir} \\
        --min-seed-length ${params.seed_finder_min_length} \\
        --tss-window ${params.seed_finder_tss_window} \\
        ${python_arg} \\
        --outdir seed_finder_work \\
        --raw-hits ${chunk_id}_raw_hits.csv \\
        --out ${chunk_id}.final.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
    END_VERSIONS
    """
}
