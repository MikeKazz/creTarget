/*
 * FLASHFRY_SCORE
 *
 * Runs FlashFry's `score` mode on discovered sites, tagging each guide
 * with on-target and off-target scoring metrics — including
 * `jostandsantos`, the specificity score developed specifically for
 * CRISPRi (as opposed to cutting Cas9), which matters for this tool's
 * use case.
 */

process FLASHFRY_SCORE {
    label 'flashfry_score'
    tag "${discover_output.simpleName}"

    container 'TODO/flashfry-image:1.15'

    input:
    path discover_output
    path flashfry_db
    val  max_mismatch
    val  scoring_metrics
    val  flashfry_jar

    output:
    path "scored.output", emit: scored_output

    script:
    """
    java -Xmx7g -jar ${flashfry_jar} \\
        score \\
        --input ${discover_output} \\
        --output scored.output \\
        --scoringMetrics ${scoring_metrics} \\
        --database ${flashfry_db} \\
        --maxMismatch ${max_mismatch}
    """
}
