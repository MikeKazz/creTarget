#!/bin/bash
#
# Verify a completed creTarget run that included --run_seed_finder.
#
#   bash verify_seed_run.sh /path/to/results_dir
#
# Read-only. Everything streams, so it is safe on a 50+ GB table: no file is
# ever loaded into memory and nothing is sorted on disk.
#
# The check that matters most is STRAND BALANCE. The scanner this release
# replaced indexed each guide under one orientation only, so it reported hits
# on a single strand and missed about half of all true sites while exiting 0.
# A one-sided strand column is the visible symptom, and it is the one thing
# you cannot detect by reading totals.

set -uo pipefail

RESULTS="${1:-}"
if [ -z "$RESULTS" ] || [ ! -d "$RESULTS" ]; then
    echo "usage: bash verify_seed_run.sh /path/to/results_dir" >&2
    exit 2
fi

RANKED="${RESULTS}/guides.ranked.tsv"
FINAL="${RESULTS}/guides.final.tsv"
HITDIR="${RESULTS}/seed_finder"

fail=0
note() { printf '  %-44s %s\n' "$1" "$2"; }
bad()  { printf '  FAIL %-39s %s\n' "$1" "$2"; fail=1; }

echo "=============================================================="
echo " creTarget seed-screening run verification"
echo " $RESULTS"
echo "=============================================================="
echo

# ---------------------------------------------------------------- files ----
echo "-- Files"
for f in "$RANKED" "$FINAL"; do
    if [ -r "$f" ]; then
        note "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    else
        bad "$(basename "$f")" "missing"
    fi
done
shards=$(find "$HITDIR" -name '*_raw_hits.csv' 2>/dev/null | wc -l | tr -d ' ')
note "seed shards" "$shards hit table(s)"
[ "$shards" -gt 0 ] || bad "seed shards" "no hit tables found in $HITDIR"
echo

# ------------------------------------------------------------ row counts ----
echo "-- Row counts"
ranked_rows=$(( $(wc -l < "$RANKED") - 1 ))
final_rows=$(( $(wc -l < "$FINAL") - 1 ))
note "guides in guides.ranked.tsv" "$(printf "%'d" "$ranked_rows")"
note "guides in guides.final.tsv"  "$(printf "%'d" "$final_rows")"
if [ "$ranked_rows" -ne "$final_rows" ]; then
    bad "row counts" "ranked and final disagree -- the seed stage lost or duplicated rows"
fi

# Exactly one header line in each: collectFile's skip/keepHeader going wrong
# shows up as an interior header, which silently corrupts any downstream parse.
rh=$(head -1 "$RANKED" | cut -f1)
for f in "$RANKED" "$FINAL"; do
    n=$(grep -c "^${rh}	" "$f")
    [ "$n" -eq 1 ] || bad "$(basename "$f") headers" "found $n header lines, expected 1"
done
note "header lines" "1 in each (ok)"
echo

# --------------------------------------------------------------- columns ----
echo "-- Columns"
missing=""
for c in seedfinder_hit_count seedfinder_closest_gene \
         seedfinder_closest_tss_distance seedfinder_closest_seed_length; do
    head -1 "$FINAL" | tr '\t' '\n' | grep -qx "$c" || missing="$missing $c"
done
if [ -n "$missing" ]; then
    bad "seed columns" "missing:$missing"
else
    note "seed columns" "all four present"
fi
echo

# ------------------------------------------------- screened vs unscreened ----
# Guides excluded by the hard filters are not screened and must be BLANK, not
# a confident 0. Mixing those up would overstate how much was checked.
echo "-- Screening coverage"
awk -F'\t' '
NR==1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "seedfinder_hit_count")   c = i
        if ($i == "passes_hard_filters")    p = i
    }
    next
}
{
    blank = ($c == "")
    if (p) {
        pass = ($p == "True" || $p == "true" || $p == "1")
        if (pass && blank)        mismatch_pass_blank++
        else if (!pass && !blank) mismatch_fail_filled++
    }
    if (blank) { unscreened++ }
    else {
        screened++
        if ($c + 0 > 0) { withhits++; hits += $c } else { nohits++ }
    }
}
END {
    printf "  %-44s %'"'"'d\n", "screened (passed hard filters)", screened
    printf "  %-44s %'"'"'d\n", "not screened (blank columns)",   unscreened
    if (screened > 0) {
        printf "  %-44s %'"'"'d  (%.1f%%)\n", "screened with >=1 seed hit", withhits, 100*withhits/screened
        printf "  %-44s %'"'"'d  (%.1f%%)\n", "screened with no seed hit",  nohits,   100*nohits/screened
        printf "  %-44s %'"'"'d\n", "total hits counted in the table", hits
    }
    if (mismatch_pass_blank)  printf "  FAIL %-39s %d rows passed the filters but were not screened\n", "coverage", mismatch_pass_blank
    if (mismatch_fail_filled) printf "  FAIL %-39s %d rows failed the filters but were screened\n", "coverage", mismatch_fail_filled
    if (mismatch_pass_blank || mismatch_fail_filled) exit 1
}' "$FINAL" || fail=1
echo

# ----------------------------------------------------------- strand balance ----
echo "-- Strand balance in the hit tables   <-- the regression check"
awk -F',' '
FNR==1 {
    s = 0
    # Tolerate CRLF: releases before 0.3.2 wrote RFC 4180 line endings, so the
    # last header field arrives as "ucsc_url\r" and a plain == misses it.
    for (i = 1; i <= NF; i++) { gsub(/\r$/, "", $i); if ($i == "strand") s = i }
    if (!s) { printf "  FAIL no strand column in %s\n", FILENAME; bad = 1 }
    next
}
s { v = $s; gsub(/\r$/, "", v); n[v]++; total++ }
END {
    if (bad) exit 1
    printf "  %-44s %'"'"'d\n", "total seed hits", total
    for (k in n) printf "  %-44s %'"'"'d  (%.1f%%)\n", "strand " k, n[k], 100*n[k]/total
    if (!(("+") in n) || !(("-") in n)) {
        print "  FAIL hits are on ONE strand only."
        print "       A scanner that indexes guides under a single orientation"
        print "       looks exactly like this and misses ~half of all true sites."
        exit 1
    }
    ratio = n["+"] / n["-"]
    if (ratio < 0.5 || ratio > 2.0)
        printf "  WARN +/- ratio is %.2f, which is more lopsided than expected\n", ratio
    else
        printf "  %-44s %.2f (balanced)\n", "+/- ratio", ratio
}' $(find "$HITDIR" -name '*_raw_hits.csv' | sort) || fail=1
echo

# ------------------------------------------------------- TSS annotation ----
# A column that is NEVER multi-valued across a genome-wide run is the
# signature of the all-TSSs-in-window behaviour regressing to nearest-only.
# Bidirectional promoters are ~10% of human genes, so some hits must sit in
# two promoters at once.
echo "-- TSS annotation"
awk -F',' '
FNR==1 {
    g = d = 0
    for (i = 1; i <= NF; i++) {
        gsub(/\r$/, "", $i)
        if ($i == "tss_genes_in_window")     g = i
        if ($i == "tss_distances_in_window") d = i
    }
    if (!g || !d) { missing = 1 }
    next
}
g {
    total++
    n = split($g, parts, ";")
    if (n > 1) multi++
    if (n > maxg) maxg = n
    if (split($d, dd, ";") != n) mismatch++
    if ($d + 0 < 0) negative++
}
END {
    if (missing) {
        print "  FAIL hit table has no tss_genes_in_window / tss_distances_in_window"
        print "       (a reference or release older than 0.5.0)"
        exit 1
    }
    printf "  %-44s %d\n", "most genes nominated by one hit", maxg
    printf "  %-44s %d  (%.2f%%)\n", "hits sitting in >1 promoter", multi+0, 100*multi/total
    printf "  %-44s %d  (%.2f%%)\n", "hits upstream of the nearest TSS", negative+0, 100*negative/total
    if (mismatch) {
        printf "  FAIL %-39s %d hit(s) have unequal gene/distance list lengths\n", "tss lists", mismatch
        exit 1
    }
    if (multi == 0) {
        print "  FAIL no hit nominates more than one gene."
        print "       Bidirectional promoters are ~10% of human genes, so a"
        print "       genome-wide run should find plenty. This looks like the"
        print "       annotation having regressed to nearest-TSS-only."
        exit 1
    }
    if (negative == 0) {
        print "  WARN no hit has a negative distance -- tss_distance should be"
        print "       SIGNED (negative upstream). An all-positive column"
        print "       suggests an unsigned value."
    }
}' $(find "$HITDIR" -name '*_raw_hits.csv' | sort) || fail=1
echo

# -------------------------------------------------------- guide_id scoping ----
# Within a shard, ids are unique by construction (the row index). Across
# shards they are unique only if each id carries its own shard's prefix --
# that is the v0.3.1 fix, and this verifies it held on real data without
# needing to sort tens of millions of ids.
echo "-- guide_id uniqueness across shards"
if [ "$shards" -le 1 ]; then
    # An unsharded run deliberately omits the prefix: the row index alone is
    # unique within one file, and short ids are easier to read. Checking for a
    # prefix here would fail every single-task run.
    note "single shard" "ids need no prefix (ok)"
else
    offenders=0
    for f in $(find "$HITDIR" -name '*_raw_hits.csv' | sort); do
        expect=$(basename "$f" _raw_hits.csv)
        wrong=$(awk -F',' -v want="$expect" '
            FNR==1 { for (i=1;i<=NF;i++) { gsub(/\r$/, "", $i); if ($i=="guide_id") g=i }; next }
            g && index($g, want "_") != 1 { bad++ }
            END { print bad + 0 }' "$f")
        if [ "$wrong" -ne 0 ]; then
            bad "$(basename "$f")" "$wrong id(s) not prefixed with '$expect'"
            offenders=$((offenders + 1))
        fi
    done
    if [ "$offenders" -eq 0 ]; then
        note "every id carries its shard prefix" "ids are globally unique (ok)"
    fi
fi
echo

# -------------------------------------------------------------- provenance ----
echo "-- Provenance"
for f in software_versions.yml run_manifest.json pipeline_info/trace.txt; do
    [ -r "${RESULTS}/${f}" ] && note "$f" "present" || note "$f" "absent"
done
SEEDREF="${CRETARGET_CACHE_DIR:-$HOME/.cretarget}/seed_reference/provenance.json"
if [ -r "$SEEDREF" ]; then
    note "seed reference provenance" "$SEEDREF"
    grep -E '"(gencode_release|tss_count|promoter_window_count)"' "$SEEDREF" \
        | sed 's/^/      /'
else
    note "seed reference provenance" "not found at $SEEDREF"
fi
echo

echo "=============================================================="
if [ "$fail" -eq 0 ]; then
    echo " ALL CHECKS PASSED"
else
    echo " SOME CHECKS FAILED -- see FAIL lines above"
fi
echo "=============================================================="
exit "$fail"
