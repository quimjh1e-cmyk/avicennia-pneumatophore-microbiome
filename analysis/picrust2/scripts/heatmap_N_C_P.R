## ============================================================
## Heatmap: N-cycle, C-cycle, P-cycle, S-cycle pathways across 24 samples
## Project: TFM — Avicennia germinans pneumatophore microbiome
## Input:   output/pathways_described.tsv
## Output:  figures/heatmap_N_C_P.pdf / .png
## ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(grid)
library(data.table)

# ── 0. Output directory ─────────────────────────────────────
dir.create("figures", showWarnings = FALSE)

# ── 1. Load pathway abundance table ─────────────────────────
# Use data.table::fread — read.table silently truncates this file at ~337 rows
# Merge two runs: nfix run supplies N2FIX-PWY/PWY-7576 (--no_gap_fill required for
# single-reaction pathways); original run supplies all other pathways with gap-fill intact.
pw_nfix <- as.data.frame(
  fread("output/pathways_nfix/pathways_described.tsv", sep = "\t", header = TRUE, check.names = FALSE),
  stringsAsFactors = FALSE)
pw_orig <- as.data.frame(
  fread("output/pathways_described.tsv", sep = "\t", header = TRUE, check.names = FALSE),
  stringsAsFactors = FALSE)
rownames(pw_nfix) <- pw_nfix[[1]]; pw_nfix <- pw_nfix[, -1]
rownames(pw_orig) <- pw_orig[[1]]; pw_orig <- pw_orig[, -1]
# Rows present only in nfix run (the two N-fixation pathways)
nfix_only <- setdiff(rownames(pw_nfix), rownames(pw_orig))
pw_dt <- rbind(pw_orig, pw_nfix[nfix_only, , drop = FALSE])
pw <- pw_dt

# Column 1 is 'description' — separate it
descr  <- pw[, 1, drop = FALSE]
abund  <- pw[, -1]           # numeric matrix: pathways × samples

# ── 2. Define focal pathways ────────────────────────────────
pathway_meta <- tribble(
  ~pathway_id,        ~label,                                        ~cycle,
  # Nitrogen
  "N2FIX-PWY",           "Nitrogen fixation I\n(ferredoxin)",         "Nitrogen",
  "DENITRIFICATION-PWY", "Nitrate reduction I\n(denitrification)",    "Nitrogen",
  "PWY490-3",            "Nitrate reduction VI\n(assimilatory)",      "Nitrogen",
  "PWY-4984",            "Urea cycle",                                "Nitrogen",
  # Carbon
  "CALVIN-PWY",          "Calvin-Benson-\nBassham cycle",             "Carbon",
  "P23-PWY",             "Reductive TCA\ncycle I",                    "Carbon",
  "PWY-5392",            "Reductive TCA\ncycle II",                   "Carbon",
  "CODH-PWY",            "Reductive acetyl-\nCoA pathway",            "Carbon",
  # Phosphorus
  "PWY-4702",            "Phytate\ndegradation I",                    "Phosphorus",
  "PWY0-1533",           "Methylphosphonate\ndegradation I",          "Phosphorus",
  "PWY-7399",            "Methylphosphonate\ndegradation II",         "Phosphorus",
  "PWY-7805",            "(Aminomethyl)phosphonate\ndegradation",     "Phosphorus",
  # Sulfur
  "SO4ASSIM-PWY",        "Assimilatory sulfate\nreduction I",         "Sulfur",
  "PWY1ZNC-1",           "Assimilatory sulfate\nreduction IV",        "Sulfur",
  "SULFATE-CYS-PWY",     "Sulfate assimilation &\ncysteine biosynthesis", "Sulfur",
  "PWY-6676",            "Sulfide oxidation\n(phototrophic)",         "Sulfur"
)

# Check which pathways are present in the table
missing <- setdiff(pathway_meta$pathway_id, rownames(abund))
if (length(missing) > 0) {
  message("WARNING – pathways not found in output: ", paste(missing, collapse=", "))
  pathway_meta <- pathway_meta[!pathway_meta$pathway_id %in% missing, ]
}

# Subset abundance to focal pathways (keep row order = pathway_meta order)
mat <- as.matrix(abund[pathway_meta$pathway_id, ])

# ── 3. Sample ordering & annotation ─────────────────────────
# Expected sample names: FCBAG1–6, FEZAG1–6, SCBAG1–6, SEZAG1–6
sample_order <- c(
  paste0("FCBAG", 1:6),
  paste0("FEZAG", 1:6),
  paste0("SCBAG", 1:6),
  paste0("SEZAG", 1:6)
)
# Only keep columns present in data (guards against name differences)
sample_order <- intersect(sample_order, colnames(mat))
# If that misses some, fall back to detected order
if (length(sample_order) < ncol(mat)) {
  message("Some sample names not matched — using original column order")
  sample_order <- colnames(mat)
}
mat <- mat[, sample_order]

# Build sample annotation (top bar)
sample_ann <- data.frame(
  Season = ifelse(grepl("^F", sample_order), "February", "September"),
  Site   = ifelse(grepl("CBA", sample_order), "Cola de Ballena", "Estero Zacatecas"),
  row.names = sample_order
)

# ── 4. Row (pathway) annotation ─────────────────────────────
row_ann <- data.frame(
  Cycle = pathway_meta$cycle,
  row.names = pathway_meta$label
)

# ── 5. z-score transformation (row-wise) ────────────────────
zscore_rows <- function(m) {
  t(apply(m, 1, function(x) {
    s <- sd(x)
    if (s == 0) rep(0, length(x)) else (x - mean(x)) / s
  }))
}
mat_z <- zscore_rows(mat)
rownames(mat_z) <- pathway_meta$label   # pretty labels for y-axis
colnames(mat_z) <- sample_order

# ── 6. Colour palettes ───────────────────────────────────────
# Annotation colours
ann_colors <- list(
  Season = c(February   = "#00BFC4",   # D3 cyan
             September  = "#C77CFF"),   # D3 purple
  Site   = c("Cola de Ballena"    = "#7CAE00",  # D3 green
             "Estero Zacatecas"   = "#F8766D"), # D3 salmon
  Cycle  = c(Nitrogen   = "#619CFF",   # D3 blue
             Carbon     = "#8C564B",   # D3 brown
             Phosphorus = "#FF7F0E",   # D3 orange
             Sulfur     = "#E377C2")   # D3 pink
)

# Diverging colour scale for z-scores
pal_breaks <- seq(-2.5, 2.5, length.out = 101)
pal_colors <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

# ── 7. Gaps between pathway groups ───────────────────────────
# Gaps between pathway groups: N=4, C=4, P=4, S=4
row_gaps <- c(4, 8, 12)    # after N, C, P blocks
# Gap after every 6 columns (groups of 6 replicates)
col_gaps <- c(6, 12, 18)

# ── 8. Draw heatmap ─────────────────────────────────────────
cairo_pdf("figures/heatmap_N_C_P_S.pdf", width = 12, height = 9, family = "sans")

pheatmap(
  mat_z,
  color             = pal_colors,
  breaks            = pal_breaks,
  cluster_rows      = FALSE,       # preserve biological grouping
  cluster_cols      = FALSE,       # preserve sample grouping
  annotation_col    = sample_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  gaps_row          = row_gaps,
  gaps_col          = col_gaps,
  fontsize          = 8,
  fontsize_row      = 7.5,
  fontsize_col      = 7.5,
  angle_col         = 45,
  border_color      = "grey90",
  cellwidth         = 22,
  cellheight        = 28,
  legend_breaks     = c(-2, -1, 0, 1, 2),
  legend_labels     = c("-2", "-1", "0", "1", "2"),
  main              = "Predicted functional pathway abundance\n(z-score per pathway)"
)

dev.off()
message("Saved: figures/heatmap_N_C_P_S.pdf")

# ── 9. PNG version (300 DPI) ─────────────────────────────────
png("figures/heatmap_N_C_P_S.png", width = 12, height = 9,
    units = "in", res = 300)

pheatmap(
  mat_z,
  color             = pal_colors,
  breaks            = pal_breaks,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = sample_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  gaps_row          = row_gaps,
  gaps_col          = col_gaps,
  fontsize          = 8,
  fontsize_row      = 7.5,
  fontsize_col      = 7.5,
  angle_col         = 45,
  border_color      = "grey90",
  cellwidth         = 22,
  cellheight        = 28,
  legend_breaks     = c(-2, -1, 0, 1, 2),
  legend_labels     = c("-2", "-1", "0", "1", "2"),
  main              = "Predicted functional pathway abundance\n(z-score per pathway)"
)

dev.off()
message("Saved: figures/heatmap_N_C_P_S.png")

# ── 10. Print raw abundance summary for reference ────────────
cat("\n── Raw abundance summary (selected pathways) ──\n")
raw_summary <- as.data.frame(mat) %>%
  rownames_to_column("pathway_id") %>%
  left_join(pathway_meta[, c("pathway_id","label","cycle")], by = "pathway_id") %>%
  mutate(
    mean_abund = rowMeans(select(., starts_with("F"), starts_with("S"))),
    .before = cycle
  ) %>%
  select(cycle, label, mean_abund, pathway_id)
print(raw_summary, row.names = FALSE)

# ============================================================
# FIGURE 2 — Nitrogenase gene cluster (KO-level)
# Nitrogen fixation pathways (N2FIX-PWY, PWY-7576) were absent
# from MinPath output; core nif KOs are shown as supplementary
# gene-level evidence consistent with Rivularia abundance.
# Input:  output/KO_described.tsv
# Output: figures/heatmap_nif_KO.pdf / .png
# ============================================================

# ── F2-1. Load KO table ──────────────────────────────────────
ko_dt <- fread("output/KO_described.tsv", sep = "\t", header = TRUE,
               check.names = FALSE)
ko <- as.data.frame(ko_dt, stringsAsFactors = FALSE)
rownames(ko) <- ko[[1]]
ko <- ko[, -1]   # drop KO id column; col 1 is now 'description'
ko_descr <- ko[, 1, drop = FALSE]
ko_abund <- ko[, -1]

# ── F2-2. Define focal nitrogenase KOs ───────────────────────
nif_meta <- tribble(
  ~ko_id,      ~label,
  "ko:K02588",  "nifH — nitrogenase\niron protein",
  "ko:K02586",  "nifD — nitrogenase MoFe\nprotein α chain",
  "ko:K02591",  "nifK — nitrogenase MoFe\nprotein β chain",
  "ko:K02587",  "nifE — MoFe cofactor\nsynthesis",
  "ko:K02592",  "nifN — MoFe cofactor\nsynthesis NifN",
  "ko:K02593",  "nifT — nitrogen\nfixation protein",
  "ko:K02595",  "nifW — nitrogenase\nstabilizing protein",
  "ko:K02596",  "nifX — nitrogen\nfixation protein"
)

missing_ko <- setdiff(nif_meta$ko_id, rownames(ko_abund))
if (length(missing_ko) > 0) {
  message("WARNING – KOs not found: ", paste(missing_ko, collapse = ", "))
  nif_meta <- nif_meta[!nif_meta$ko_id %in% missing_ko, ]
}

# ── F2-3. Build matrix, apply same sample order ───────────────
mat_ko <- as.matrix(ko_abund[nif_meta$ko_id, sample_order])

# ── F2-4. z-score (row-wise) ─────────────────────────────────
mat_ko_z <- zscore_rows(mat_ko)
rownames(mat_ko_z) <- nif_meta$label
colnames(mat_ko_z) <- sample_order

# ── F2-5. Draw heatmap (same styling as Figure 1) ────────────
cairo_pdf("figures/heatmap_nif_KO.pdf", width = 12, height = 6, family = "sans")

pheatmap(
  mat_ko_z,
  color             = pal_colors,
  breaks            = pal_breaks,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = sample_ann,
  annotation_colors = ann_colors[c("Season", "Site")],
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  gaps_col          = col_gaps,
  fontsize          = 8,
  fontsize_row      = 7.5,
  fontsize_col      = 7.5,
  angle_col         = 45,
  border_color      = "grey90",
  cellwidth         = 22,
  cellheight        = 28,
  legend_breaks     = c(-2, -1, 0, 1, 2),
  legend_labels     = c("-2", "-1", "0", "1", "2"),
  main              = "Predicted nitrogenase gene cluster abundance\n(z-score per KO; pathway-level inference not available — see Methods)"
)

dev.off()
message("Saved: figures/heatmap_nif_KO.pdf")

png("figures/heatmap_nif_KO.png", width = 12, height = 6,
    units = "in", res = 300)

pheatmap(
  mat_ko_z,
  color             = pal_colors,
  breaks            = pal_breaks,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = sample_ann,
  annotation_colors = ann_colors[c("Season", "Site")],
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  gaps_col          = col_gaps,
  fontsize          = 8,
  fontsize_row      = 7.5,
  fontsize_col      = 7.5,
  angle_col         = 45,
  border_color      = "grey90",
  cellwidth         = 22,
  cellheight        = 28,
  legend_breaks     = c(-2, -1, 0, 1, 2),
  legend_labels     = c("-2", "-1", "0", "1", "2"),
  main              = "Predicted nitrogenase gene cluster abundance\n(z-score per KO; pathway-level inference not available — see Methods)"
)

dev.off()
message("Saved: figures/heatmap_nif_KO.png")

# ── F2-6. Print group means for reporting ────────────────────
cat("\n── nif KO mean abundance by group ──\n")
nif_summary <- as.data.frame(mat_ko) %>%
  rownames_to_column("ko_id") %>%
  left_join(nif_meta, by = "ko_id") %>%
  mutate(
    FCBA = rowMeans(select(., starts_with("FCBA"))),
    FEZA = rowMeans(select(., starts_with("FEZA"))),
    SCBA = rowMeans(select(., starts_with("SCBA"))),
    SEZA = rowMeans(select(., starts_with("SEZA")))
  ) %>%
  select(ko_id, label, FCBA, FEZA, SCBA, SEZA)
print(nif_summary, row.names = FALSE)
