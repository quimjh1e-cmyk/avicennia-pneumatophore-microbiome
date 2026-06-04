# ============================================================
# ggpicrust2 — Visualization & KEGG Analysis
# Project : TFM — Avicennia germinans pneumatophore microbiome
# ============================================================
#
# CONTEXT
# -------
# Tasks 1 & 2 are COMPLETE in microeco (Picrust_R_analysis.R):
#   - microtables built: func_pw, func_ko, func_ec
#   - 8-method DA run for season / site / season_site
#   - consensus synthesis CSVs saved in Output/4.PICRUSt2/
#
# THIS SCRIPT adds ggpicrust2-specific outputs:
#   A) KEGG pathway analysis — NEW (ko2kegg_abundance; microeco can't do this)
#   B) Errorbar plots — reuse existing microeco DESeq2 results, no re-testing
#   C) PCA + heatmap visualizations for all three comparisons
#   D) Targeted N/P cycle plots
#
# Reference: Yang et al. 2023, Bioinformatics 39(8), btad470
# ============================================================


# ============================================================
# 0.  Packages
# ============================================================
# install.packages("ggpicrust2")
# install.packages(c("tidyverse", "ggprism", "patchwork"))
# BiocManager::install(c("ALDEx2", "DESeq2", "edgeR", "limma"))

library(ggpicrust2)
library(tidyverse)
library(ggplot2)
library(ggprism)
library(patchwork)
library(GGally)    # required internally by pathway_errorbar()

# Set working directory to analysis/microeco/ before running this script


# ============================================================
# 1.  Load existing microeco results
# ============================================================
# func_pw / func_ko / func_ec  : microtable objects (abundance + metadata)
# pw_desc / ko_desc / ec_desc  : named vectors of descriptions
load("Output/4.PICRUSt2/picrust2_microtables.RData")

# Convenience: extract abundance matrices and metadata
pw_abun   <- as.data.frame(func_pw$otu_table)    # 539 pathways × 24 samples
ko_abun   <- as.data.frame(func_ko$otu_table)    # 8,822 KOs    × 24 samples
ec_abun   <- as.data.frame(func_ec$otu_table)    # 2,944 ECs    × 24 samples
meta      <- func_pw$sample_table                 # 24 samples

# season_site interaction column (if not already present)
if (!"season_site" %in% colnames(meta)) {
  meta$season_site <- paste(meta$season, meta$site, sep = "_")
}

output_dir <- "Output/4.PICRUSt2"


# ============================================================
# 2.  Helper — reformat microeco DESeq2 CSV → ggpicrust2 format
# ============================================================
# ggpicrust2's pathway_errorbar() expects a daa_results data.frame with:
#   feature, p_values, p_adjust, log2FoldChange, group1, group2,
#   method, description
#
# microeco DESeq2 output columns:
#   Comparison, Taxa (= "ID|ID"), baseMean, log2FoldChange, P.unadj, P.adj, Method

# Generic loader — handles both DESeq2 (log2FoldChange) and
# metagenomeSeq (logFC) microeco output formats.
load_microeco_da <- function(csv_path, desc_vec = pw_desc) {
  df <- read.csv(csv_path, row.names = 1, check.names = FALSE)

  # Extract clean pathway/feature ID (strip duplicate "ID|ID" → "ID")
  df$feature <- sub("\\|.*", "", df$Taxa)

  # Normalise fold-change column name: metagenomeSeq uses 'logFC',
  # DESeq2 uses 'log2FoldChange'.  Both are log2-scaled.
  if ("logFC" %in% colnames(df) && !"log2FoldChange" %in% colnames(df)) {
    df <- df %>% rename(log2FoldChange = logFC)
  }

  # Parse comparison string e.g. "February - September"
  # For season_site there are multiple pairwise Comparison values; we store
  # the per-row group1/group2 so enriched_in is correct for every row.
  df <- df %>%
    mutate(
      group1 = trimws(sub(" - .*", "", Comparison)),
      group2 = trimws(sub(".* - ", "", Comparison))
    )

  df <- df %>%
    rename(
      p_values = P.unadj,
      p_adjust = P.adj,
      method   = Method
    ) %>%
    mutate(
      description = desc_vec[feature],
      # Positive log2FC = enriched in group1 (numerator of the comparison)
      enriched_in = ifelse(log2FoldChange > 0, group1, group2)
    ) %>%
    select(feature, group1, group2, enriched_in,
           p_values, p_adjust, log2FoldChange, method, description) %>%
    filter(!is.na(feature))

  df
}


# ============================================================
# 3.  Helper — load consensus to filter to robust features
# ============================================================
# Returns pathway IDs found by ≥ threshold methods
load_consensus_ids <- function(consensus_csv, min_methods = 6) {
  df <- read.csv(consensus_csv, check.names = FALSE)
  df$Pathway_ID[df$Total_Methods_Sig >= min_methods]
}


# ============================================================
# 4.  Helper — errorbar plot from reformatted results
# ============================================================
make_errorbar <- function(
    abun,           # abundance data.frame (features × samples)
    daa_df,         # reformatted daa results (ggpicrust2 format)
    meta_vec,       # named character vector: sample → group label
    pathway_type,   # "MetaCyc" | "EC" | "KEGG"
    ko_to_kegg  = FALSE,
    consensus_ids = NULL,   # if not NULL, restrict to these IDs
    p_threshold = 0.05,
    desc_vec    = NULL,     # named vector ID → description for NA backfill
    label = "plot"
) {
  # Filter to consensus-significant features only
  if (!is.null(consensus_ids)) {
    daa_df <- daa_df %>% filter(feature %in% consensus_ids)
  }

  sig <- daa_df %>% filter(p_adjust < p_threshold)
  n_sig <- nrow(sig)
  cat(label, "— significant features:", n_sig, "\n")
  if (n_sig == 0) { message("Nothing to plot for ", label); return(invisible(NULL)) }

  # Annotate: add pathway descriptions.
  # When desc_vec is supplied (MetaCyc, EC) we bypass pathway_annotation()
  # entirely — our desc_vec comes directly from PICRUSt2's *_described.tsv
  # and covers every feature ID, so there is nothing to gain from the
  # internal lookup, and it avoids the "missing value where TRUE/FALSE needed"
  # error that arises from NAs in ggpicrust2's annotation table.
  if (!is.null(desc_vec)) {
    annotated <- sig %>%
      mutate(description = {
        d <- desc_vec[feature]
        ifelse(is.na(d) | d == "", feature, d)
      })
  } else {
    annotated <- tryCatch(
      pathway_annotation(sig, pathway_type, ko_to_kegg),
      error = function(e) { message("annotation failed: ", e$message); sig }
    )
    if (!"description" %in% colnames(annotated)) annotated$description <- annotated$feature
    na_d <- is.na(annotated$description) | annotated$description == ""
    annotated$description[na_d] <- annotated$feature[na_d]
  }

  # Carry enriched_in through annotation (it may have been dropped by
  # pathway_annotation() if it ran).  Recompute from log2FoldChange.
  if (!"enriched_in" %in% colnames(annotated)) {
    annotated <- annotated %>%
      mutate(enriched_in = ifelse(log2FoldChange > 0, group1, group2))
  }

  # Build and print a direction table for the console
  direction_label <- paste0("Positive log2FC = enriched in ",
                            unique(sig$group1)[1],
                            "  (vs ", unique(sig$group2)[1], ")")
  cat("  Direction: ", direction_label, "\n", sep = "")

  top_summary <- annotated %>%
    filter(p_adjust < p_threshold) %>%
    arrange(p_adjust) %>%
    slice_head(n = 20) %>%
    select(feature, description, log2FoldChange, p_adjust, enriched_in)
  cat("  Top significant pathways:\n")
  print(as.data.frame(top_summary), row.names = FALSE)

  # Build Group factor with explicit level order so that pathway_errorbar()'s
  # internal log2(level1 / level2) matches DESeq2's convention.
  # DESeq2 encodes comparison as "group1 - group2", so group1 is the numerator
  # (positive log2FC = enriched in group1).  pathway_errorbar() computes
  # log2(level1 / level2) where level1 is factor level 1.
  # By making DESeq2's group1 = factor level 1, the two sign conventions agree.
  grp1 <- unique(annotated$group1)[1]
  grp2 <- unique(annotated$group2)[1]
  Group <- factor(meta_vec[colnames(abun)], levels = c(grp1, grp2))

  # Cap display at 30 features (pathway_errorbar max)
  MAX_FEATURES <- 30
  if (n_sig > MAX_FEATURES)
    cat("  Showing top", MAX_FEATURES, "of", n_sig, "by p_adjust\n")

  top_ids <- annotated %>%
    filter(p_adjust < p_threshold) %>%
    arrange(p_adjust) %>%
    slice_head(n = MAX_FEATURES) %>%
    pull(feature)

  p <- tryCatch(
    pathway_errorbar(
      abundance          = abun,
      daa_results        = annotated,
      Group              = Group,
      ko_to_kegg         = ko_to_kegg,
      p_values_threshold = p_threshold,
      order              = if (pathway_type == "KEGG") "pathway_class" else "group",
      select             = top_ids,
      p_value_bar        = TRUE,
      x_lab              = if (pathway_type == "KEGG") "pathway_name" else "description"
    ),
    error = function(e) { message("errorbar failed: ", e$message); NULL }
  )

  if (!is.null(p)) {
    # Patch a subtitle onto the ggplot object so the PDF is self-explanatory
    p <- p + labs(subtitle = direction_label)

    n_shown <- length(top_ids)
    height  <- max(6, n_shown * 0.45 + 3)
    ggsave(file.path(output_dir, paste0("errorbar_", label, ".pdf")),
           plot = p, width = 14, height = height)
    cat("  Saved: errorbar_", label, ".pdf\n", sep = "")

    # Save enrichment direction table alongside the plot.
    # Wrapped in tryCatch: on Windows the file may be locked (open in Excel).
    csv_path <- file.path(output_dir, paste0("enrichment_direction_", label, ".csv"))
    tryCatch(
      write.csv(
        annotated %>% filter(p_adjust < p_threshold) %>% arrange(p_adjust),
        csv_path, row.names = FALSE
      ),
      error = function(e) message(
        "  Could not save ", basename(csv_path),
        " (file may be open in another program): ", e$message
      )
    )
  }
  invisible(p)
}


# ============================================================
# 5.  Helper — PCA and heatmap
# ============================================================
make_pca_heatmap <- function(abun, meta_df, group_col, label,
                             desc_vec = NULL) {
  meta_sub <- meta_df %>%
    select(all_of(group_col)) %>%
    rename(group = 1)

  # PCA
  p_pca <- tryCatch(
    pathway_pca(abundance = abun, metadata = meta_sub, group = "group"),
    error = function(e) { message("PCA failed: ", e$message); NULL }
  )
  if (!is.null(p_pca)) {
    ggsave(file.path(output_dir, paste0("pca_", label, ".pdf")),
           p_pca, width = 7, height = 6)
    cat("  Saved: pca_", label, ".pdf\n", sep = "")
  }

  # Heatmap — top 50 most variable features
  top50 <- names(sort(apply(abun, 1, var), decreasing = TRUE)[1:50])
  abun_top <- abun[top50, ]

  # Replace row names with human-readable descriptions when available.
  # make.unique() handles the rare case of two IDs sharing a description.
  if (!is.null(desc_vec)) {
    new_names           <- desc_vec[top50]
    new_names[is.na(new_names) | new_names == ""] <- top50[is.na(new_names) | new_names == ""]
    rownames(abun_top)  <- make.unique(new_names, sep = " #")
  }

  p_hm <- tryCatch(
    pathway_heatmap(abundance = abun_top, metadata = meta_sub, group = "group"),
    error = function(e) { message("heatmap failed: ", e$message); NULL }
  )
  if (!is.null(p_hm)) {
    ggsave(file.path(output_dir, paste0("heatmap_", label, ".pdf")),
           p_hm, width = 10, height = 13)
    cat("  Saved: heatmap_", label, ".pdf\n", sep = "")
  }
  invisible(list(pca = p_pca, heatmap = p_hm))
}


# ============================================================
# SECTION A — KEGG Pathway Analysis (NEW — microeco can't do this)
# ============================================================
# ko2kegg_abundance() maps 8,822 KO IDs to KEGG pathways and
# sums their abundances. This produces a KEGG-pathway × sample table
# with KEGG class hierarchy — unavailable in microeco.

cat("\n=== A. KEGG Pathway Analysis (new) ===\n")

# ko2kegg_abundance expects: data.frame with first col named 'function'
# containing KO IDs WITHOUT "ko:" prefix, remaining cols = samples
ko_for_kegg <- ko_abun %>%
  rownames_to_column("function") %>%
  mutate(`function` = sub("^ko:", "", `function`))

# ko2kegg_abundance() only accepts a file path string, not a data.frame.
# Write to a temp TSV and pass the path.
tmp_ko <- tempfile(fileext = ".tsv")
write.table(ko_for_kegg, tmp_ko, sep = "\t", row.names = FALSE, quote = FALSE)

kegg_abun <- tryCatch(
  ko2kegg_abundance(tmp_ko),
  error = function(e) { message("ko2kegg_abundance failed: ", e$message); NULL }
)
unlink(tmp_ko)   # clean up

if (!is.null(kegg_abun)) {
  cat("KEGG pathways after conversion:", nrow(kegg_abun), "\n")

  # Save the KEGG abundance table
  write.csv(kegg_abun,
            file.path(output_dir, "KEGG_pathway_abundance.csv"))

  # Group vectors
  season_vec      <- setNames(meta$season,      rownames(meta))
  site_vec        <- setNames(meta$site,        rownames(meta))
  season_site_vec <- setNames(meta$season_site, rownames(meta))

  # Helper: annotate KEGG DA results.
  # pathway_annotation() consistently fails on our pre-converted IDs due to
  # NAs in its internal lookup table, so we skip it and build the required
  # columns directly. The feature IDs (e.g. "ko00010") are already meaningful
  # KEGG codes and serve as pathway names until a working annotation is available.
  annotate_kegg_safe <- function(daa_df) {
    ann <- daa_df
    if (!"pathway_name"  %in% colnames(ann)) ann$pathway_name  <- ann$feature
    if (!"pathway_class" %in% colnames(ann)) ann$pathway_class <- "KEGG"
    if (!"description"   %in% colnames(ann)) ann$description   <- ann$feature
    # Clean up any residual NAs
    for (col in c("pathway_name", "pathway_class", "description")) {
      bad <- is.na(ann[[col]]) | ann[[col]] == ""
      ann[[col]][bad] <- ann$feature[bad]
    }
    ann
  }

  kegg_errorbar <- function(ann_df, Group_vec, file_label) {
    n_sig <- sum(ann_df$p_adjust < 0.05, na.rm = TRUE)
    cat("  Significant KEGG pathways:", n_sig, "\n")
    if (n_sig == 0) return(invisible(NULL))

    top_ids <- ann_df %>%
      filter(p_adjust < 0.05) %>% arrange(p_adjust) %>%
      slice_head(n = 30) %>% pull(feature)

    # For KEGG (LinDA): make the factor order explicit so bar colours are
    # consistent.  LinDA uses the alphabetically first level as the reference
    # (negative coefficient = enriched in reference), so level1 = ref group.
    kegg_lvls   <- sort(unique(Group_vec))   # alphabetical = LinDA reference first
    kegg_Group  <- factor(Group_vec[colnames(kegg_abun)], levels = kegg_lvls)

    p_err <- tryCatch(
      pathway_errorbar(
        abundance          = kegg_abun,
        daa_results        = ann_df,
        Group              = kegg_Group,
        ko_to_kegg         = FALSE,
        p_values_threshold = 0.05,
        select             = top_ids,
        order              = "pathway_class",
        x_lab              = "pathway_name"
      ),
      error = function(e) { message("  KEGG errorbar failed: ", e$message); NULL }
    )
    if (!is.null(p_err))
      ggsave(file.path(output_dir, paste0("errorbar_KEGG_", file_label, ".pdf")),
             p_err, width = 14, height = max(6, length(top_ids) * 0.45 + 3))
    invisible(p_err)
  }

  # --- A1. Season ---
  cat("\nA1. KEGG × season\n")
  daa_kegg_season <- tryCatch(
    pathway_daa(abundance = kegg_abun, metadata = meta %>% select(season),
                group = "season", daa_method = "LinDA", p_adjust_method = "BH"),
    error = function(e) { message(e$message); NULL }
  )
  ann_kegg_season <- NULL
  if (!is.null(daa_kegg_season)) {
    ann_kegg_season <- annotate_kegg_safe(daa_kegg_season)
    write.csv(ann_kegg_season,
              file.path(output_dir, "DA_KEGG_season_LinDA.csv"), row.names = FALSE)
    kegg_errorbar(ann_kegg_season, season_vec, "season")
  }

  # --- A2. Site ---
  cat("\nA2. KEGG × site\n")
  daa_kegg_site <- tryCatch(
    pathway_daa(abundance = kegg_abun, metadata = meta %>% select(site),
                group = "site", daa_method = "LinDA", p_adjust_method = "BH"),
    error = function(e) { message(e$message); NULL }
  )
  ann_kegg_site <- NULL
  if (!is.null(daa_kegg_site)) {
    ann_kegg_site <- annotate_kegg_safe(daa_kegg_site)
    write.csv(ann_kegg_site,
              file.path(output_dir, "DA_KEGG_site_LinDA.csv"), row.names = FALSE)
    kegg_errorbar(ann_kegg_site, site_vec, "site")
  }

  # --- A3. PCA + heatmap ---
  # Build KEGG pathway-name lookup for heatmap labels.
  # Use ann_kegg_season if available; fall back to feature IDs.
  kegg_desc <- if (!is.null(ann_kegg_season) &&
                    "pathway_name" %in% colnames(ann_kegg_season)) {
    nm <- ann_kegg_season$pathway_name
    nm[is.na(nm)] <- ann_kegg_season$feature[is.na(nm)]
    setNames(nm, ann_kegg_season$feature)
  } else NULL

  cat("\nA3. KEGG PCA / heatmap × season_site\n")
  make_pca_heatmap(kegg_abun, meta, "season_site", "KEGG_season_site", desc_vec = kegg_desc)
  make_pca_heatmap(kegg_abun, meta, "season",      "KEGG_season",      desc_vec = kegg_desc)
  make_pca_heatmap(kegg_abun, meta, "site",        "KEGG_site",        desc_vec = kegg_desc)
}


# ============================================================
# SECTION B — MetaCyc Errorbar Plots (reuse existing DA)
# ============================================================
# Primary method for DIRECTION: DESeq2
#   - Best power/FDR in the 8-method benchmark for this dataset
#   - Raw-abundance verification confirmed DESeq2/LinDA call correct direction
#   - metagenomeSeq CSS normalization systematically reversed directions on this
#     dataset despite having better power/FDR in the benchmark
#
# Significance filter: consensus (≥6/8 methods) — metagenomeSeq's high
# sensitivity still contributes through the consensus, without trusting
# its direction calls.

cat("\n\n=== B. MetaCyc Errorbar Plots (metagenomeSeq) ===\n")

season_vec      <- setNames(meta$season,      rownames(meta))
site_vec        <- setNames(meta$site,        rownames(meta))
season_site_vec <- setNames(meta$season_site, rownames(meta))

# B1. Season
cat("\nB1. MetaCyc × season\n")
daa_pw_season  <- load_microeco_da("Output/4.PICRUSt2/DA_season_DESeq2.csv", pw_desc)
cons_pw_season <- load_consensus_ids("Output/4.PICRUSt2/DA_season_Consensus_Synthesis.csv")
make_errorbar(pw_abun, daa_pw_season, season_vec, "MetaCyc",
              consensus_ids = cons_pw_season, desc_vec = pw_desc,
              label = "MetaCyc_season")

# B2. Site
cat("\nB2. MetaCyc × site\n")
daa_pw_site  <- load_microeco_da("Output/4.PICRUSt2/DA_site_DESeq2.csv", pw_desc)
cons_pw_site <- load_consensus_ids("Output/4.PICRUSt2/DA_site_Consensus_Synthesis.csv")
make_errorbar(pw_abun, daa_pw_site, site_vec, "MetaCyc",
              consensus_ids = cons_pw_site, desc_vec = pw_desc,
              label = "MetaCyc_site")

# B3. PCA + heatmap (all three groupings)
cat("\nB3. MetaCyc PCA / heatmap\n")
make_pca_heatmap(pw_abun, meta, "season_site", "MetaCyc_season_site", desc_vec = pw_desc)
make_pca_heatmap(pw_abun, meta, "season",      "MetaCyc_season",      desc_vec = pw_desc)
make_pca_heatmap(pw_abun, meta, "site",        "MetaCyc_site",        desc_vec = pw_desc)


# ============================================================
# SECTION C — EC Number Errorbar Plots (reuse existing DA)
# ============================================================
cat("\n\n=== C. EC Number Errorbar Plots ===\n")

# C1. Season
cat("\nC1. EC × season\n")
# Note: microeco DA results for EC are in the same output folder
ec_season_path <- "Output/4.PICRUSt2/DA_season_DESeq2.csv"   # replace with EC-specific if available
# If a dedicated EC DA CSV exists, load that; otherwise note that
# EC errorbar requires EC-specific DA results. Uncomment the block
# below once EC DA CSVs are confirmed:
#
# daa_ec_season  <- load_microeco_da("Output/4.PICRUSt2/DA_EC_season_DESeq2.csv", ec_desc)
# cons_ec_season <- load_consensus_ids  ("Output/4.PICRUSt2/DA_EC_season_Consensus_Synthesis.csv")
# make_errorbar(ec_abun, daa_ec_season, season_vec, "EC",
#               consensus_ids = cons_ec_season, label = "EC_season")

# PCA + heatmap always works without DA results
cat("C2. EC PCA / heatmap\n")
make_pca_heatmap(ec_abun, meta, "season_site", "EC_season_site", desc_vec = ec_desc)
make_pca_heatmap(ec_abun, meta, "season",      "EC_season",      desc_vec = ec_desc)
make_pca_heatmap(ec_abun, meta, "site",        "EC_site",        desc_vec = ec_desc)


# ============================================================
# SECTION D — Targeted N & P Cycle Plots
# ============================================================
cat("\n\n=== D. Targeted N & P Cycle Pathway Plots ===\n")

# MetaCyc pathway IDs related to N and P cycling
nitrogen_ids <- names(pw_desc)[grepl(
  "nitrif|denitrif|nitrogen|nitrate|nitrite|ammonia|ammon|urea|N-fix|dinitrogen|glutamine.synth",
  pw_desc, ignore.case = TRUE)]

phosphorus_ids <- names(pw_desc)[grepl(
  "phosph|polyphosphate|phytate",
  pw_desc, ignore.case = TRUE)]

cat("Nitrogen-related pathways found:", length(nitrogen_ids), "\n")
cat("Phosphorus-related pathways found:", length(phosphorus_ids), "\n")

# Long-format data for boxplots
make_targeted_boxplot <- function(ids, cycle_name, abun, meta_df) {
  present <- intersect(ids, rownames(abun))
  if (length(present) == 0) {
    message(cycle_name, ": no matching pathways in abundance table")
    return(invisible(NULL))
  }

  df <- abun[present, ] %>%
    rownames_to_column("Pathway_ID") %>%
    pivot_longer(-Pathway_ID, names_to = "sample", values_to = "abundance") %>%
    left_join(meta_df %>% rownames_to_column("sample"), by = "sample") %>%
    mutate(
      description  = pw_desc[Pathway_ID],
      # Use description as facet label; fall back to ID if missing
      facet_label  = ifelse(is.na(description) | description == "",
                            Pathway_ID,
                            str_wrap(description, width = 40))
    )

  p <- ggplot(df, aes(x = season, y = abundance, fill = season)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
    facet_grid(facet_label ~ site, scales = "free_y") +
    scale_fill_manual(values = c("February" = "seagreen", "September" = "#F28E2B")) +
    labs(
      title = paste(cycle_name, "cycle — predicted pathway abundance"),
      x     = NULL, y = "Predicted abundance", fill = "Season"
    ) +
    theme_prism(base_size = 9) +
    theme(
      strip.text.y = element_text(size = 6, angle = 0, hjust = 0),
      axis.text.x  = element_text(angle = 30, hjust = 1)
    )

  h <- max(8, length(present) * 1.8)
  ggsave(file.path(output_dir,
                   paste0("targeted_", tolower(cycle_name), "_pathways.pdf")),
         p, width = 10, height = h, limitsize = FALSE)
  cat("  Saved:", paste0("targeted_", tolower(cycle_name), "_pathways.pdf"), "\n")
  invisible(p)
}

make_targeted_boxplot(nitrogen_ids,  "Nitrogen",   pw_abun, meta)
make_targeted_boxplot(phosphorus_ids, "Phosphorus", pw_abun, meta)


# ============================================================
# SECTION E — Summary table: consensus-significant pathways
# ============================================================
cat("\n\n=== E. Consensus summary tables ===\n")

summarise_consensus <- function(csv_path, desc_vec, label) {
  df <- read.csv(csv_path, check.names = FALSE) %>%
    mutate(description = desc_vec[Pathway_ID]) %>%
    arrange(desc(Total_Methods_Sig))

  write.csv(df, file.path(output_dir,
                           paste0("consensus_annotated_", label, ".csv")),
            row.names = FALSE)
  cat(" ", label, "— total sig:", nrow(df[df$Total_Methods_Sig >= 1, ]),
      "| consensus (≥6):", nrow(df[df$Total_Methods_Sig >= 6, ]), "\n")
  invisible(df)
}

summarise_consensus("Output/4.PICRUSt2/DA_season_Consensus_Synthesis.csv",
                    pw_desc, "season")
summarise_consensus("Output/4.PICRUSt2/DA_site_Consensus_Synthesis.csv",
                    pw_desc, "site")
summarise_consensus("Output/4.PICRUSt2/DA_season_site_Consensus_Synthesis.csv",
                    pw_desc, "season_site")


# ============================================================
# Done
# ============================================================
cat("\n\nAll outputs written to:", output_dir, "\n")
cat("Files created:\n")
cat(paste0("  ", list.files(output_dir, pattern = "\\.(pdf|csv)$")), sep = "\n")
