process ANNOTATE {
    tag "${chunk_id}"

    input:
    tuple val(chunk_id), path(scored)
    path region_metadata
    path atac_narrowpeak

    output:
    tuple val(chunk_id), path("${chunk_id}.annotated.tsv"), emit: annotated
    path "versions.yml", emit: versions

    script:
    def atac_arg = atac_narrowpeak.name != 'NO_ATAC_FILE' ? "--atac-narrowpeak ${atac_narrowpeak}" : ''
    """
    cretarget annotate \\
        --scored ${scored} \\
        --region-metadata ${region_metadata} \\
        --blast-radius-bp ${params.blast_radius_bp} \\
        ${atac_arg} \\
        --out ${chunk_id}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
    END_VERSIONS
    """
}
