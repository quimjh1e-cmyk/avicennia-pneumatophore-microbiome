#!/usr/bin/env bash
#
# PICRUSt2 pipeline for 16S amplicon data
#
# Predicts metagenomic functional potential from ASV data using PICRUSt2.
# Input: ASV abundance table (.tsv) and representative sequences (.fasta)
# Output: Predicted MetaCyc pathway, KEGG Ortholog (KO), and EC abundances
#
# Usage:
#   ./run_picrust2_pipeline.sh -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [-n NSTI_CUTOFF]
#
# Requirements:
#   - conda environment with picrust2 and biom-format installed
#   - Input files: asv_table.tsv (tab-separated, first column = ASV IDs)
#                  rep_seqs.fasta (representative sequences matching ASV IDs)
#
# Example:
#   conda activate picrust2
#   ./run_picrust2_pipeline.sh -i ./input -o ./output -t 4 -n 2

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
THREADS=4
NSTI_CUTOFF=2
INPUT_DIR=""
OUTPUT_DIR=""

# ── Parse arguments ───────────────────────────────────────────────────────────
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [-n NSTI_CUTOFF]"
    echo ""
    echo "  -i  Directory containing asv_table.tsv and rep_seqs.fasta"
    echo "  -o  Output directory for PICRUSt2 results"
    echo "  -t  Number of threads (default: 4)"
    echo "  -n  NSTI cutoff for filtering unreliable predictions (default: 2)"
    exit 1
}

while getopts "i:o:t:n:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        n) NSTI_CUTOFF="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

ASV_TABLE="$INPUT_DIR/asv_table.tsv"
REP_SEQS="$INPUT_DIR/rep_seqs.fasta"
BIOM_TABLE="$INPUT_DIR/asv_table.biom"

for f in "$ASV_TABLE" "$REP_SEQS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/picrust2_run.log"

echo "========================================"
echo "PICRUSt2 Pipeline"
echo "========================================"
echo "Input:   $INPUT_DIR"
echo "Output:  $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "NSTI cutoff: $NSTI_CUTOFF"
echo "Log:     $LOG_FILE"
echo "========================================"

# ── Step 1: Convert ASV table to BIOM format ─────────────────────────────────
echo ""
echo "[Step 1/5] Converting ASV table to BIOM format..."

# Ensure header uses "#OTU ID" (BIOM requirement)
HEADER=$(head -1 "$ASV_TABLE")
if [[ "$HEADER" != "#OTU ID"* ]]; then
    echo "  Fixing header: replacing first column name with '#OTU ID'"
    TMP_TABLE="$INPUT_DIR/asv_table_biom_ready.tsv"
    sed '1s/^[^\t]*/\#OTU ID/' "$ASV_TABLE" > "$TMP_TABLE"
else
    TMP_TABLE="$ASV_TABLE"
fi

biom convert \
    -i "$TMP_TABLE" \
    -o "$BIOM_TABLE" \
    --table-type="OTU table" \
    --to-hdf5

SAMPLE_COUNT=$(biom summarize-table -i "$BIOM_TABLE" 2>/dev/null | grep "Num samples" | awk '{print $NF}')
ASV_COUNT=$(biom summarize-table -i "$BIOM_TABLE" 2>/dev/null | grep "Num observations" | awk '{print $NF}')
echo "  BIOM table: $ASV_COUNT ASVs x $SAMPLE_COUNT samples"

# Clean up temp file if created
[[ "$TMP_TABLE" != "$ASV_TABLE" ]] && rm -f "$TMP_TABLE"

# ── Step 2: Run PICRUSt2 full pipeline ────────────────────────────────────────
echo ""
echo "[Step 2/5] Running PICRUSt2 pipeline (this may take 15-45 min)..."

picrust2_pipeline.py \
    -s "$REP_SEQS" \
    -i "$BIOM_TABLE" \
    -o "$OUTPUT_DIR" \
    -p "$THREADS" \
    --stratified \
    --per_sequence_contrib \
    --verbose \
    2>&1 | tee "$LOG_FILE"

# ── Step 3: NSTI quality control ─────────────────────────────────────────────
echo ""
echo "[Step 3/5] NSTI quality control..."

NSTI_FILE="$OUTPUT_DIR/marker_predicted_and_nsti.tsv.gz"
if [[ -f "$NSTI_FILE" ]]; then
    NSTI_STATS=$(zcat "$NSTI_FILE" | awk -F'\t' '
        NR>1 {
            sum+=$NF; count++;
            if($NF>max) max=$NF;
            if(min=="" || $NF<min) min=$NF;
            if($NF>cutoff) excluded++
        }
        END {
            printf "  ASVs placed:    %d\n  Mean NSTI:      %.4f\n  Min NSTI:       %.4f\n  Max NSTI:       %.4f\n  ASVs > %.1f:     %d\n",
                count, sum/count, min, max, cutoff, excluded+0
        }' cutoff="$NSTI_CUTOFF")
    echo "$NSTI_STATS"
    echo "$NSTI_STATS" > "$OUTPUT_DIR/nsti_summary.txt"
else
    echo "  WARNING: NSTI file not found at $NSTI_FILE"
fi

# ── Step 4: Add human-readable descriptions ──────────────────────────────────
echo ""
echo "[Step 4/5] Adding pathway/KO/EC descriptions..."

add_descriptions.py \
    -i "$OUTPUT_DIR/pathways_out/path_abun_unstrat.tsv.gz" \
    -m METACYC \
    -o "$OUTPUT_DIR/pathways_described.tsv"
echo "  MetaCyc pathways: $(tail -n +2 "$OUTPUT_DIR/pathways_described.tsv" | wc -l) features"

add_descriptions.py \
    -i "$OUTPUT_DIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    -m KO \
    -o "$OUTPUT_DIR/KO_described.tsv"
echo "  KO orthologs:     $(tail -n +2 "$OUTPUT_DIR/KO_described.tsv" | wc -l) features"

add_descriptions.py \
    -i "$OUTPUT_DIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    -m EC \
    -o "$OUTPUT_DIR/EC_described.tsv"
echo "  EC enzymes:       $(tail -n +2 "$OUTPUT_DIR/EC_described.tsv" | wc -l) features"

# ── Step 5: Summary ──────────────────────────────────────────────────────────
echo ""
echo "[Step 5/5] Pipeline complete."
echo "========================================"
echo "Key output files:"
echo "  $OUTPUT_DIR/pathways_described.tsv  — MetaCyc pathway abundances"
echo "  $OUTPUT_DIR/KO_described.tsv        — KEGG Ortholog abundances"
echo "  $OUTPUT_DIR/EC_described.tsv        — Enzyme Commission abundances"
echo "  $OUTPUT_DIR/nsti_summary.txt        — NSTI quality summary"
echo "  $OUTPUT_DIR/picrust2_run.log        — Full run log"
echo ""
echo "Next: import the described TSVs into R for downstream analysis."
echo "========================================"
