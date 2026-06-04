#!/usr/bin/env bash
#
# QIIME2 pipeline for 16S amplicon analysis (paired-end, DADA2, SILVA)
#
# Designed for reproducible 16S rRNA amplicon workflows. All study-specific
# parameters are configurable via command-line flags or a config file.
# Defaults correspond to the A. germinans pneumatophore study (TFM UOC).
#
# Usage:
#   conda activate qiime2-amplicon-2026.1
#   cd /path/to/project
#   bash scripts/qiime2/qiime2_pipeline.sh [OPTIONS]
#
# Options:
#   --config FILE        Source a config file (bash variables) instead of flags
#   --primer-f SEQ       Forward primer sequence        [CCTACGGGNGGCWGCAG]
#   --primer-r SEQ       Reverse primer sequence        [GACTACHVGGGTATCTAATCC]
#   --trunc-f INT        DADA2 forward truncation       [279]
#   --trunc-r INT        DADA2 reverse truncation       [200]
#   --maxee-f NUM        DADA2 max expected errors fwd  [2]
#   --maxee-r NUM        DADA2 max expected errors rev  [3]
#   --rarefaction INT    Rarefaction sampling depth      [21000]
#   --silva-ver STR      SILVA version                   [138.2]
#   --threads INT        CPU threads                     [4]
#   --raw-dir DIR        Directory with FASTQ files      [data/raw]
#   --qiime-dir DIR      Output directory for artifacts  [data/qiime2]
#   --metadata FILE      Sample metadata TSV             [data/metadata/metadata.tsv]
#   --fwd-suffix STR     Forward read filename suffix    [_1.fastq.gz]
#   --rev-suffix STR     Reverse read filename suffix    [_2.fastq.gz]
#   --help               Show this help message
#
# Config file example (pipeline.conf):
#   PRIMER_F="GTGYCAGCMGCCGCGGTAA"
#   PRIMER_R="CCGYCAATTYMTTTRAGTTT"
#   TRUNC_F=220
#   TRUNC_R=160
#   RAREFACTION_DEPTH=10000
#
# Each step checks whether its output already exists and skips if so.
# To force a rerun, delete the relevant output file(s).

set -euo pipefail

# ============================================================
# DEFAULT PARAMETERS
# ============================================================
PRIMER_F="CCTACGGGNGGCWGCAG"
PRIMER_R="GACTACHVGGGTATCTAATCC"
TRUNC_F=279
TRUNC_R=200
MAXEE_F=2
MAXEE_R=3
RAREFACTION_DEPTH=21000
SILVA_VER="138.2"
THREADS=4
RAW_DIR="data/raw"
QIIME_DIR="data/qiime2"
META_DIR="data/metadata"
METADATA="data/metadata/metadata.tsv"
FWD_SUFFIX="_1.fastq.gz"
REV_SUFFIX="_2.fastq.gz"

# ============================================================
# PARSE ARGUMENTS
# ============================================================
show_help() {
    sed -n '/^# Usage:/,/^[^#]/p' "$0" | head -n -1 | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)       source "$2"; shift 2 ;;
        --primer-f)     PRIMER_F="$2"; shift 2 ;;
        --primer-r)     PRIMER_R="$2"; shift 2 ;;
        --trunc-f)      TRUNC_F="$2"; shift 2 ;;
        --trunc-r)      TRUNC_R="$2"; shift 2 ;;
        --maxee-f)      MAXEE_F="$2"; shift 2 ;;
        --maxee-r)      MAXEE_R="$2"; shift 2 ;;
        --rarefaction)  RAREFACTION_DEPTH="$2"; shift 2 ;;
        --silva-ver)    SILVA_VER="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --raw-dir)      RAW_DIR="$2"; shift 2 ;;
        --qiime-dir)    QIIME_DIR="$2"; shift 2 ;;
        --metadata)     METADATA="$2"; shift 2 ;;
        --fwd-suffix)   FWD_SUFFIX="$2"; shift 2 ;;
        --rev-suffix)   REV_SUFFIX="$2"; shift 2 ;;
        --help|-h)      show_help ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ============================================================
# VALIDATION
# ============================================================
if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: metadata file not found: $METADATA" >&2
    echo "       Run this script from the project root directory." >&2
    exit 1
fi

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: raw data directory not found: $RAW_DIR" >&2
    exit 1
fi

n_samples=$(find "$RAW_DIR" -name "*${FWD_SUFFIX}" | wc -l)
if [[ "$n_samples" -eq 0 ]]; then
    echo "ERROR: no forward reads (*${FWD_SUFFIX}) found in $RAW_DIR" >&2
    exit 1
fi

echo "=== QIIME2 16S amplicon pipeline ==="
echo "  Primers:     $PRIMER_F / $PRIMER_R"
echo "  DADA2:       trunc $TRUNC_F/$TRUNC_R, maxEE $MAXEE_F/$MAXEE_R"
echo "  Rarefaction: $RAREFACTION_DEPTH"
echo "  SILVA:       $SILVA_VER"
echo "  Samples:     $n_samples (from $RAW_DIR)"
echo "  Threads:     $THREADS"
echo ""

# ============================================================
# GENERATE MANIFEST
# ============================================================
# QIIME2 manifest requires absolute paths; generate at runtime from $PWD.
MANIFEST_DIR="$(dirname "$METADATA")"
MANIFEST="$MANIFEST_DIR/manifest_abs.tsv"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"
for fwd in "$RAW_DIR"/*"$FWD_SUFFIX"; do
    sample=$(basename "$fwd" "$FWD_SUFFIX")
    rev="$RAW_DIR/${sample}${REV_SUFFIX}"
    if [[ ! -f "$rev" ]]; then
        echo "ERROR: reverse read not found for $sample: $rev" >&2
        exit 1
    fi
    echo -e "${sample}\t$(realpath "$fwd")\t$(realpath "$rev")"
done >> "$MANIFEST"

mkdir -p "$QIIME_DIR"

# ============================================================
# STEP 1: Import demultiplexed paired-end reads
# ============================================================
echo "=== Step 1: Import ==="
if [[ ! -f "$QIIME_DIR/demux-paired-end.qza" ]]; then
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-format PairedEndFastqManifestPhred33V2 \
        --input-path "$MANIFEST" \
        --output-path "$QIIME_DIR/demux-paired-end.qza"
else
    echo "  Skipping (demux-paired-end.qza exists)"
fi

# ============================================================
# STEP 2: Pre-trimming quality summary
# ============================================================
echo "=== Step 2: Pre-trimming quality summary ==="
if [[ ! -f "$QIIME_DIR/demux-paired-end-summary.qzv" ]]; then
    qiime demux summarize \
        --i-data "$QIIME_DIR/demux-paired-end.qza" \
        --o-visualization "$QIIME_DIR/demux-paired-end-summary.qzv"
else
    echo "  Skipping (demux-paired-end-summary.qzv exists)"
fi

# ============================================================
# STEP 3: Primer trimming with Cutadapt
# ============================================================
# Cutadapt performs sequence-based primer matching (handles degenerate bases).
# --p-discard-untrimmed removes reads where primers are not found.
# --p-no-indels requires exact-length primer matches.
echo "=== Step 3: Cutadapt primer trimming ==="
if [[ ! -f "$QIIME_DIR/trimmed-seqs.qza" ]]; then
    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "$QIIME_DIR/demux-paired-end.qza" \
        --p-front-f "$PRIMER_F" \
        --p-front-r "$PRIMER_R" \
        --p-discard-untrimmed \
        --p-no-indels \
        --p-cores "$THREADS" \
        --o-trimmed-sequences "$QIIME_DIR/trimmed-seqs.qza" \
        --verbose
else
    echo "  Skipping (trimmed-seqs.qza exists)"
fi

# ============================================================
# STEP 4: Post-trimming quality summary
# ============================================================
echo "=== Step 4: Post-trimming quality summary ==="
if [[ ! -f "$QIIME_DIR/trimmed-seqs-summary.qzv" ]]; then
    qiime demux summarize \
        --i-data "$QIIME_DIR/trimmed-seqs.qza" \
        --o-visualization "$QIIME_DIR/trimmed-seqs-summary.qzv"
else
    echo "  Skipping (trimmed-seqs-summary.qzv exists)"
fi

# ============================================================
# STEP 5: DADA2 denoising
# ============================================================
# Truncation must satisfy: trunc_f + trunc_r >= amplicon_length + min_overlap (12 bp).
# Adjust --trunc-f and --trunc-r based on your quality profiles (Step 2/4 summaries).
echo "=== Step 5: DADA2 denoising ==="
if [[ ! -f "$QIIME_DIR/feature-table.qza" ]]; then
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$QIIME_DIR/trimmed-seqs.qza" \
        --p-trunc-len-f "$TRUNC_F" \
        --p-trunc-len-r "$TRUNC_R" \
        --p-max-ee-f "$MAXEE_F" \
        --p-max-ee-r "$MAXEE_R" \
        --p-n-threads "$THREADS" \
        --o-table "$QIIME_DIR/feature-table.qza" \
        --o-representative-sequences "$QIIME_DIR/rep-seqs.qza" \
        --o-denoising-stats "$QIIME_DIR/dada2-stats.qza" \
        --verbose
    # --o-base-transition-stats is available in QIIME2 v2026.1 but omitted
    # here; it outputs DADA2's learned error model (base transition rates)
    # which is useful for diagnostics but not required for downstream analysis.
else
    echo "  Skipping (feature-table.qza exists)"
fi

# ============================================================
# STEP 6: DADA2 summary
# ============================================================
echo "=== Step 6: Feature table summary ==="
if [[ ! -f "$QIIME_DIR/feature-table-summary.qzv" ]]; then
    qiime feature-table summarize \
        --i-table "$QIIME_DIR/feature-table.qza" \
        --m-metadata-file "$METADATA" \
        --o-feature-frequencies "$QIIME_DIR/feature-frequencies.qza" \
        --o-sample-frequencies "$QIIME_DIR/sample-frequencies.qza" \
        --o-summary "$QIIME_DIR/feature-table-summary.qzv"
else
    echo "  Skipping (feature-table-summary.qzv exists)"
fi

# ============================================================
# STEP 7: Train SILVA classifier on target region
# ============================================================
# Uses RESCRIPt to download SILVA, extract the primer region, and train
# a Naive Bayes classifier de novo.
SILVA_PREFIX="silva-${SILVA_VER}"

echo "=== Step 7a: Download SILVA ${SILVA_VER} ==="
if [[ ! -f "$QIIME_DIR/${SILVA_PREFIX}-seqs.qza" ]]; then
    qiime rescript get-silva-data \
        --p-version "$SILVA_VER" \
        --p-target 'SSURef_NR99' \
        --o-silva-sequences "$QIIME_DIR/${SILVA_PREFIX}-seqs.qza" \
        --o-silva-taxonomy "$QIIME_DIR/${SILVA_PREFIX}-tax.qza"
else
    echo "  Skipping (${SILVA_PREFIX}-seqs.qza exists)"
fi

echo "=== Step 7b: Reverse-transcribe RNA to DNA ==="
if [[ ! -f "$QIIME_DIR/${SILVA_PREFIX}-dna-seqs.qza" ]]; then
    qiime rescript reverse-transcribe \
        --i-rna-sequences "$QIIME_DIR/${SILVA_PREFIX}-seqs.qza" \
        --o-dna-sequences "$QIIME_DIR/${SILVA_PREFIX}-dna-seqs.qza"
else
    echo "  Skipping (${SILVA_PREFIX}-dna-seqs.qza exists)"
fi

echo "=== Step 7c: Extract primer region ==="
if [[ ! -f "$QIIME_DIR/${SILVA_PREFIX}-ref-seqs.qza" ]]; then
    qiime feature-classifier extract-reads \
        --i-sequences "$QIIME_DIR/${SILVA_PREFIX}-dna-seqs.qza" \
        --p-f-primer "$PRIMER_F" \
        --p-r-primer "$PRIMER_R" \
        --p-min-length 200 \
        --p-max-length 500 \
        --p-n-jobs "$THREADS" \
        --o-reads "$QIIME_DIR/${SILVA_PREFIX}-ref-seqs.qza"
else
    echo "  Skipping (${SILVA_PREFIX}-ref-seqs.qza exists)"
fi

echo "=== Step 7d: Train Naive Bayes classifier ==="
if [[ ! -f "$QIIME_DIR/${SILVA_PREFIX}-classifier.qza" ]]; then
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "$QIIME_DIR/${SILVA_PREFIX}-ref-seqs.qza" \
        --i-reference-taxonomy "$QIIME_DIR/${SILVA_PREFIX}-tax.qza" \
        --o-classifier "$QIIME_DIR/${SILVA_PREFIX}-classifier.qza"
else
    echo "  Skipping (${SILVA_PREFIX}-classifier.qza exists)"
fi

# ============================================================
# STEP 8: Classify ASVs
# ============================================================
echo "=== Step 8: Taxonomic classification ==="
if [[ ! -f "$QIIME_DIR/taxonomy.qza" ]]; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "$QIIME_DIR/${SILVA_PREFIX}-classifier.qza" \
        --i-reads "$QIIME_DIR/rep-seqs.qza" \
        --p-n-jobs "$THREADS" \
        --o-classification "$QIIME_DIR/taxonomy.qza"
else
    echo "  Skipping (taxonomy.qza exists)"
fi

echo "=== Step 8b: Taxonomy visualization ==="
if [[ ! -f "$QIIME_DIR/taxonomy-summary.qzv" ]]; then
    qiime metadata tabulate \
        --m-input-file "$QIIME_DIR/taxonomy.qza" \
        --o-visualization "$QIIME_DIR/taxonomy-summary.qzv"
else
    echo "  Skipping (taxonomy-summary.qzv exists)"
fi

echo "=== Step 8c: Taxonomy barplot (pre-filtering) ==="
if [[ ! -f "$QIIME_DIR/taxa-barplot.qzv" ]]; then
    qiime taxa barplot \
        --i-table "$QIIME_DIR/feature-table.qza" \
        --i-taxonomy "$QIIME_DIR/taxonomy.qza" \
        --m-metadata-file "$METADATA" \
        --o-visualization "$QIIME_DIR/taxa-barplot.qzv"
else
    echo "  Skipping (taxa-barplot.qzv exists)"
fi

# ============================================================
# STEP 9: Remove chloroplast and mitochondrial reads
# ============================================================
# Universal 16S primers co-amplify chloroplast and mitochondrial rRNA
# from host tissue and eukaryotic epiphytes. These are removed before
# downstream prokaryotic community analysis.
echo "=== Step 9: Remove chloroplast and mitochondrial reads ==="
if [[ ! -f "$QIIME_DIR/table-no-chloro-mito.qza" ]]; then
    qiime taxa filter-table \
        --i-table "$QIIME_DIR/feature-table.qza" \
        --i-taxonomy "$QIIME_DIR/taxonomy.qza" \
        --p-exclude "Chloroplast,Mitochondria" \
        --o-filtered-table "$QIIME_DIR/table-no-chloro-mito.qza"
else
    echo "  Skipping (table-no-chloro-mito.qza exists)"
fi

if [[ ! -f "$QIIME_DIR/rep-seqs-no-chloro-mito.qza" ]]; then
    qiime taxa filter-seqs \
        --i-sequences "$QIIME_DIR/rep-seqs.qza" \
        --i-taxonomy "$QIIME_DIR/taxonomy.qza" \
        --p-exclude "Chloroplast,Mitochondria" \
        --o-filtered-sequences "$QIIME_DIR/rep-seqs-no-chloro-mito.qza"
else
    echo "  Skipping (rep-seqs-no-chloro-mito.qza exists)"
fi

echo "=== Step 9b: Filtered barplot and summary ==="
if [[ ! -f "$QIIME_DIR/taxa-barplot-filtered.qzv" ]]; then
    qiime taxa barplot \
        --i-table "$QIIME_DIR/table-no-chloro-mito.qza" \
        --i-taxonomy "$QIIME_DIR/taxonomy.qza" \
        --m-metadata-file "$METADATA" \
        --o-visualization "$QIIME_DIR/taxa-barplot-filtered.qzv"
else
    echo "  Skipping (taxa-barplot-filtered.qzv exists)"
fi

if [[ ! -f "$QIIME_DIR/table-no-chloro-mito-summary.qzv" ]]; then
    qiime feature-table summarize \
        --i-table "$QIIME_DIR/table-no-chloro-mito.qza" \
        --m-metadata-file "$METADATA" \
        --o-feature-frequencies "$QIIME_DIR/feature-frequencies-filtered.qza" \
        --o-sample-frequencies "$QIIME_DIR/sample-frequencies-filtered.qza" \
        --o-summary "$QIIME_DIR/table-no-chloro-mito-summary.qzv"
else
    echo "  Skipping (table-no-chloro-mito-summary.qzv exists)"
fi

# ============================================================
# STEP 10: Phylogenetic tree construction
# ============================================================
# Required for phylogeny-based diversity metrics (UniFrac).
# MAFFT alignment -> mask variable positions -> FastTree -> midpoint root.
echo "=== Step 10a: MAFFT alignment ==="
if [[ ! -f "$QIIME_DIR/aligned-rep-seqs.qza" ]]; then
    qiime alignment mafft \
        --i-sequences "$QIIME_DIR/rep-seqs-no-chloro-mito.qza" \
        --p-n-threads "$THREADS" \
        --o-alignment "$QIIME_DIR/aligned-rep-seqs.qza"
else
    echo "  Skipping (aligned-rep-seqs.qza exists)"
fi

echo "=== Step 10b: Mask highly variable positions ==="
if [[ ! -f "$QIIME_DIR/masked-aligned-rep-seqs.qza" ]]; then
    qiime alignment mask \
        --i-alignment "$QIIME_DIR/aligned-rep-seqs.qza" \
        --o-masked-alignment "$QIIME_DIR/masked-aligned-rep-seqs.qza"
else
    echo "  Skipping (masked-aligned-rep-seqs.qza exists)"
fi

echo "=== Step 10c: FastTree ==="
if [[ ! -f "$QIIME_DIR/unrooted-tree.qza" ]]; then
    qiime phylogeny fasttree \
        --i-alignment "$QIIME_DIR/masked-aligned-rep-seqs.qza" \
        --p-n-threads "$THREADS" \
        --o-tree "$QIIME_DIR/unrooted-tree.qza"
else
    echo "  Skipping (unrooted-tree.qza exists)"
fi

echo "=== Step 10d: Midpoint rooting ==="
if [[ ! -f "$QIIME_DIR/rooted-tree.qza" ]]; then
    qiime phylogeny midpoint-root \
        --i-tree "$QIIME_DIR/unrooted-tree.qza" \
        --o-rooted-tree "$QIIME_DIR/rooted-tree.qza"
else
    echo "  Skipping (rooted-tree.qza exists)"
fi

# ============================================================
# STEP 11: Core diversity analysis
# ============================================================
# Rarefaction depth should be set to the minimum post-filtering sample depth
# (check Step 9b summary). Produces alpha diversity, beta diversity
# (Bray-Curtis, Jaccard, weighted/unweighted UniFrac), and PCoA ordinations.
echo "=== Step 11: Core diversity metrics ==="
if [[ ! -d "$QIIME_DIR/core-metrics-results" ]]; then
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny "$QIIME_DIR/rooted-tree.qza" \
        --i-table "$QIIME_DIR/table-no-chloro-mito.qza" \
        --p-sampling-depth "$RAREFACTION_DEPTH" \
        --m-metadata-file "$METADATA" \
        --p-n-jobs-or-threads "$THREADS" \
        --output-dir "$QIIME_DIR/core-metrics-results"
else
    echo "  Skipping (core-metrics-results/ exists)"
fi

echo ""
echo "=== Pipeline complete ==="
echo "Key outputs:"
echo "  Feature table (filtered): $QIIME_DIR/table-no-chloro-mito.qza"
echo "  Rep sequences (filtered): $QIIME_DIR/rep-seqs-no-chloro-mito.qza"
echo "  Taxonomy:                 $QIIME_DIR/taxonomy.qza"
echo "  Rooted tree:              $QIIME_DIR/rooted-tree.qza"
echo "  Core diversity:           $QIIME_DIR/core-metrics-results/"
echo ""
echo "Next: export artifacts for downstream analysis"
