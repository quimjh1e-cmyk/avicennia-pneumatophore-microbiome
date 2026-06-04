#!/usr/bin/env bash
# Search for a taxon (family or genus) across priority DA methods for SEASON comparison
# Usage: ./search_da_season.sh <taxon_name>
# Example: ./search_da_season.sh Prolixibacteraceae
#          ./search_da_season.sh Rivularia

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DA_DIR="$SCRIPT_DIR/../analysis/microeco/Output/1.Amplicon/Stage6_Diff_abund"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <taxon_name>"
    echo "  Searches family or genus name across metagenomeSeq, DESeq2, ANCOM-BC2, and LEfSe (season)"
    exit 1
fi

TAXON="$1"
SEP="$(printf '=%.0s' {1..100})"

echo "$SEP"
echo "SEASON DA SEARCH: $TAXON"
echo "Sign conventions: metagenomeSeq logFC negative = Feb | DESeq2 log2FC negative = Feb | ANCOM-BC2 lfc negative = Feb"
echo "$SEP"

# --- metagenomeSeq (ASV level) ---
echo ""
echo ">>> metagenomeSeq (Comparison: Feb - Sep | negative logFC = Feb-enriched, positive = Sep-enriched)"
echo "    Columns: Row, Comparison, Taxa, logFC, se, P.unadj, P.adj, Significance, Method"
FILE="$DA_DIR/Diff_abund_ASV_season_metagenomeSeq.csv"
MATCHES=$(grep -i "$TAXON" "$FILE" || true)
if [[ -n "$MATCHES" ]]; then
    SIG=$(echo "$MATCHES" | grep -c -E '"\*+"|"\*"' || true)
    TOTAL=$(echo "$MATCHES" | wc -l)
    echo "    Found: $TOTAL ASVs ($SIG significant)"
    echo "$MATCHES" | while IFS= read -r line; do
        sig_field=$(echo "$line" | awk -F',' '{print $(NF-1)}')
        if echo "$sig_field" | grep -q '\*'; then
            echo "  * $line"
        else
            echo "    $line"
        fi
    done
else
    echo "    No matches."
fi

# --- DESeq2 (ASV level) ---
echo ""
echo ">>> DESeq2 (Comparison: Sep - Feb | positive log2FC = Sep-enriched, negative = Feb-enriched)"
echo "    Columns: Row, Comparison, Taxa, baseMean, log2FoldChange, lfcSE, stat, P.unadj, P.adj, Significance, Method"
FILE="$DA_DIR/Diff_abund_ASV_season_DESeq2.csv"
MATCHES=$(grep -i "$TAXON" "$FILE" || true)
if [[ -n "$MATCHES" ]]; then
    SIG=$(echo "$MATCHES" | grep -c -E '"\*+"|"\*"' || true)
    TOTAL=$(echo "$MATCHES" | wc -l)
    echo "    Found: $TOTAL ASVs ($SIG significant)"
    echo "$MATCHES" | while IFS= read -r line; do
        sig_field=$(echo "$line" | awk -F',' '{print $(NF-1)}')
        if echo "$sig_field" | grep -q '\*'; then
            echo "  * $line"
        else
            echo "    $line"
        fi
    done
else
    echo "    No matches."
fi

# --- ANCOM-BC2 (ASV level, filter to season factor only) ---
echo ""
echo ">>> ANCOM-BC2 (Factor: seasonSeptember | positive lfc = Sep-enriched, negative = Feb-enriched)"
echo "    Columns: Row, Taxa, Factors, lfc, se, W, p, P.adj, diff, passed_ss, Significance, Method"
FILE="$DA_DIR/Diff_abund_ASV_season_ancombc2.csv"
MATCHES=$(grep -i "$TAXON" "$FILE" | grep '"seasonSeptember"' || true)
if [[ -n "$MATCHES" ]]; then
    SIG=$(echo "$MATCHES" | grep -c -E '"\*+"|"\*"' || true)
    TOTAL=$(echo "$MATCHES" | wc -l)
    echo "    Found: $TOTAL ASVs ($SIG significant)"
    echo "$MATCHES" | while IFS= read -r line; do
        sig_field=$(echo "$line" | awk -F',' '{print $(NF-1)}')
        if echo "$sig_field" | grep -q '\*'; then
            echo "  * $line"
        else
            echo "    $line"
        fi
    done
else
    echo "    No matches."
fi

# --- LEfSe Family level ---
echo ""
echo ">>> LEfSe - Family level (Group = enriched season)"
echo "    Columns: Row, Comparison, Taxa, Method, Group, LDA, P.unadj, P.adj, Significance"
FILE="$DA_DIR/Diff_abund_Family_season_lefse.csv"
MATCHES=$(grep -i "$TAXON" "$FILE" || true)
if [[ -n "$MATCHES" ]]; then
    echo "$MATCHES"
else
    echo "    No matches."
fi

# --- LEfSe allTax level (includes genus) ---
echo ""
echo ">>> LEfSe - All taxonomy levels (Group = enriched season)"
echo "    Columns: Row, Comparison, Taxa, Method, Group, LDA, P.unadj, P.adj, Significance"
FILE="$DA_DIR/Diff_abund_allTax_season_lefse.csv"
MATCHES=$(grep -i "$TAXON" "$FILE" || true)
if [[ -n "$MATCHES" ]]; then
    echo "$MATCHES"
else
    echo "    No matches."
fi

echo ""
echo "$SEP"
echo "SUMMARY"
echo "$SEP"

# LEfSe summary
echo ""
echo "LEfSe enrichment:"
for lefse_file in \
    "Family|$DA_DIR/Diff_abund_Family_season_lefse.csv" \
    "allTax|$DA_DIR/Diff_abund_allTax_season_lefse.csv"; do

    IFS='|' read -r level file <<< "$lefse_file"
    matches=$(grep -i "$TAXON" "$file" 2>/dev/null || true)
    if [[ -n "$matches" ]]; then
        while IFS= read -r line; do
            group=$(echo "$line" | awk -F',' '{gsub(/"/, "", $5); print $5}')
            lda=$(echo "$line" | awk -F',' '{gsub(/"/, "", $6); printf "%.2f", $6}')
            sig=$(echo "$line" | awk -F',' '{gsub(/"/, "", $9); print $9}')
            taxa=$(echo "$line" | awk -F',' '{gsub(/"/, "", $3); n=split($3, a, "|"); print a[n]}')
            echo "  LEfSe $level: $taxa → $group (LDA=$lda $sig)"
        done <<< "$matches"
    else
        echo "  LEfSe $level: no match"
    fi
done

# ASV-level direction tally
echo ""
echo "Direction tally (significant ASVs only):"
for method_file in \
    "metagenomeSeq|$DA_DIR/Diff_abund_ASV_season_metagenomeSeq.csv|4|neg=Feb,pos=Sep" \
    "DESeq2|$DA_DIR/Diff_abund_ASV_season_DESeq2.csv|5|pos=Sep,neg=Feb" \
    "ANCOM-BC2|$DA_DIR/Diff_abund_ASV_season_ancombc2.csv|4|pos=Sep,neg=Feb"; do

    IFS='|' read -r method file col convention <<< "$method_file"

    if [[ "$method" == "ANCOM-BC2" ]]; then
        sig_lines=$(grep -i "$TAXON" "$file" 2>/dev/null | grep '"seasonSeptember"' | grep -E '"\*+"|"\*"' || true)
    else
        sig_lines=$(grep -i "$TAXON" "$file" 2>/dev/null | grep -E '"\*+"|"\*"' || true)
    fi

    if [[ -n "$sig_lines" ]]; then
        pos=0; neg=0
        while IFS= read -r line; do
            val=$(echo "$line" | awk -F',' -v c="$col" '{gsub(/"/, "", $c); print $c}')
            if awk "BEGIN{exit ($val >= 0) ? 0 : 1}" 2>/dev/null; then
                ((pos++)) || true
            else
                ((neg++)) || true
            fi
        done <<< "$sig_lines"
        echo "  $method ($convention): ${pos} positive, ${neg} negative"
    else
        echo "  $method: no significant ASVs"
    fi
done
