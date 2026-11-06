process RANK {
    tag "${chunk_id}"
    label 'process_single'

    input:
    tuple val(chunk_id), path(annotated)

    output:
    tuple val(chunk_id), path("${chunk_id}.ranked.tsv"), emit: ranked
    path "versions.yml", emit: versions

    script:
    """
    cretarget rank \\
        --annotated ${annotated} \\
        --max-mismatch ${params.max_mismatch} \\
        --out ${chunk_id}.ranked.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
    END_VERSIONS
    """
}
