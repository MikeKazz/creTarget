process FLASHFRY_DISCOVER_SCORE {
    tag "${chunk.baseName}"

    input:
    path chunk
    path flashfry_jar
    tuple path(flashfry_db), path(flashfry_db_header)

    output:
    tuple val("${chunk.baseName}"), path("${chunk.baseName}.scored.tsv"), emit: scored
    path "versions.yml", emit: versions

    script:
    /*
     * Heap sizes come from the actual allocation rather than a constant, so
     * a retry at 2x memory gets a 2x heap for free.
     *
     * discover holds the off-target index largely OFF the Java heap, so it
     * gets the smaller fraction of the allocation despite being the heavier
     * step -- the remainder is headroom the JVM needs for the memory-mapped
     * index plus OS overhead. Overshooting here is the classic way to get
     * an OOM-kill on a node that looked like it had plenty of memory.
     */
    def total_gb     = task.memory ? task.memory.toGiga() : 8
    def discover_gb  = Math.max(2, (total_gb * 0.5) as int)
    def score_gb     = Math.max(2, (total_gb * 0.6) as int)
    """
    java -Xmx${discover_gb}g -jar ${flashfry_jar} \\
        discover \\
        --database ${flashfry_db} \\
        --fasta ${chunk} \\
        --output ${chunk.baseName}.discover.tsv \\
        --maxMismatch ${params.max_mismatch}

    java -Xmx${score_gb}g -jar ${flashfry_jar} \\
        score \\
        --input ${chunk.baseName}.discover.tsv \\
        --output ${chunk.baseName}.scored.tsv \\
        --scoringMetrics ${params.scoring_metrics} \\
        --database ${flashfry_db} \\
        --maxMismatch ${params.max_mismatch}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n1 | sed 's/.*version "\\(.*\\)".*/\\1/')
        flashfry: \$(basename ${flashfry_jar})
    END_VERSIONS
    """
}
