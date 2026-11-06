process PREP_REGIONS {
    tag "prep"

    input:
    path snp_bed
    path enhancer_bed
    tuple path(fasta), path(fai)

    output:
    path "regions.fasta"        , emit: fasta
    path "regions.metadata.tsv" , emit: metadata
    path "versions.yml"         , emit: versions

    script:
    def snp_arg   = snp_bed.name      != 'NO_SNP_FILE' ? "--snp-bed ${snp_bed}"           : ''
    def enh_arg   = enhancer_bed.name != 'NO_ENH_FILE' ? "--enhancer-bed ${enhancer_bed}" : ''
    def dedup_arg = params.dedup_blast_radius ? '--dedup-blast-radius' : ''
    """
    cretarget prep-regions \\
        ${snp_arg} \\
        ${enh_arg} \\
        --fasta ${fasta} \\
        --flank-bp ${params.flank_bp} \\
        --snp-window-bp ${params.snp_window_bp} \\
        --blast-radius-bp ${params.blast_radius_bp} \\
        ${dedup_arg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cretarget: \$(cretarget --version | sed 's/creTarget //')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
