#!/usr/bin/env Rscript
#
# pathway_taxa.R
# Map PICRUSt2 MetaCyc pathway predictions to contributing taxa.
#
# Given a list of pathway IDs, this script:
#   1. Resolves pathway -> reactions -> EC numbers (via PICRUSt2 mapping files)
#   2. Queries the EC-level pred_metagenome_contrib.tsv.gz for ASV contributions
#   3. Joins with QIIME2 taxonomy
#   4. Outputs ranked contribution tables at ASV, Genus, and Family levels
#   5. Optionally aggregates by sample groups (from metadata)
#
# Dependencies: data.table (only)
#
# Usage:
#   Rscript pathway_taxa.R \
#     --pathways P23-PWY,CODH-PWY,DENITRIFICATION-PWY \
#     --ec-contrib /path/to/EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
#     --taxonomy /path/to/taxonomy.tsv \
#     --output-dir results/ \
#     [--metadata /path/to/metadata.tsv --group-col treatment] \
#     [--samples FCBAG1,FCBAG2,FCBAG3]  \
#     [--path2rxn /path/to/metacyc_path2rxn_struc_filt_pro.txt] \
#     [--rxn2ec /path/to/metacyc_rxn_to_level4ec.tsv]
#
# The taxonomy file should be a QIIME2-exported taxonomy.tsv with columns:
#   Feature ID <tab> Taxon [<tab> Confidence]
# where Taxon is semicolon-delimited with rank prefixes, e.g.:
#   d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...
# Any prefix style works (d__/k__/D_0__/p__/etc.) — the script strips them.

suppressPackageStartupMessages(library(data.table))

# ── Argument parsing (base R, no optparse needed) ────────────────────────────

parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  args <- list()
  i <- 1L
  while (i <= length(raw)) {
    if (startsWith(raw[i], "--")) {
      key <- sub("^--", "", raw[i])
      if (i < length(raw) && !startsWith(raw[i + 1L], "--")) {
        args[[key]] <- raw[i + 1L]
        i <- i + 2L
      } else {
        args[[key]] <- TRUE
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  args
}

args <- parse_args()

# Validate required arguments
always_required <- c("ec-contrib", "taxonomy", "output-dir")
has_pathways <- !is.null(args[["pathways"]])
has_ecs      <- !is.null(args[["ecs"]])
missing <- setdiff(always_required, names(args))

if (length(missing) > 0L || (!has_pathways && !has_ecs)) {
  cat("pathway_taxa.R — Map PICRUSt2 pathways or ECs to contributing taxa\n\n")
  cat("Required arguments:\n")
  cat("  --pathways      Comma-separated MetaCyc pathway IDs\n")
  cat("        OR\n")
  cat("  --ecs           Comma-separated EC numbers (e.g. 1.18.6.1,1.2.7.4)\n")
  cat("                  Use this for pathways missing from the mapping file,\n")
  cat("                  or to query specific enzymes directly.\n")
  cat("  --ec-contrib    Path to EC_metagenome_out/pred_metagenome_contrib.tsv.gz\n")
  cat("  --taxonomy      Path to QIIME2 taxonomy.tsv\n")
  cat("  --output-dir    Output directory\n\n")
  cat("Optional arguments:\n")
  cat("  --metadata      Sample metadata TSV (first column = sample ID)\n")
  cat("  --group-col     Column in metadata for group comparisons\n")
  cat("  --samples       Comma-separated sample IDs to restrict analysis\n")
  cat("  --path2rxn      Override: pathway-to-reaction mapping file\n")
  cat("  --rxn2ec        Override: reaction-to-EC mapping file\n")
  cat("  --top-n         Number of top taxa to report per level (default: 20)\n")
  cat("  --label         Label for --ecs output files (default: 'custom_ECs')\n\n")
  cat("Examples:\n")
  cat("  # By pathway\n")
  cat("  Rscript pathway_taxa.R --pathways P23-PWY,CODH-PWY ...\n\n")
  cat("  # By EC (e.g. nitrogenase = EC:1.18.6.1)\n")
  cat("  Rscript pathway_taxa.R --ecs 1.18.6.1 --label nitrogenase ...\n\n")
  if (length(missing) > 0L)
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  stop("Provide either --pathways or --ecs.", call. = FALSE)
}

# Parse input mode
ec_contrib_path <- args[["ec-contrib"]]
tax_path     <- args[["taxonomy"]]
out_dir      <- args[["output-dir"]]
meta_path    <- args[["metadata"]]
group_col    <- args[["group-col"]]
sample_ids   <- if (!is.null(args[["samples"]])) trimws(unlist(strsplit(args[["samples"]], ","))) else NULL
top_n        <- as.integer(ifelse(is.null(args[["top-n"]]), 20L, args[["top-n"]]))

# Two modes: --pathways (resolve via mapping files) or --ecs (direct query)
direct_ec_mode <- has_ecs && !has_pathways
if (has_pathways && has_ecs) {
  cat("  NOTE: Both --pathways and --ecs provided. Using --pathways (--ecs ignored).\n")
  direct_ec_mode <- FALSE
}

if (direct_ec_mode) {
  direct_ecs <- trimws(unlist(strsplit(args[["ecs"]], ",")))
  direct_ecs <- sub("^EC:", "", direct_ecs, ignore.case = TRUE)
  ec_label   <- ifelse(is.null(args[["label"]]), "custom_ECs", args[["label"]])
} else {
  pathway_ids <- trimws(unlist(strsplit(args[["pathways"]], ",")))
}

# ── Step 1: Resolve target EC numbers ────────────────────────────────────────

if (direct_ec_mode) {
  # Direct EC mode: skip pathway resolution entirely
  cat(sprintf("\n[1/5] Direct EC mode: %d EC numbers provided\n", length(direct_ecs)))
  cat(sprintf("  ECs: %s\n", paste(direct_ecs, collapse = ", ")))
  cat(sprintf("  Output label: %s\n", ec_label))
  # Build pw_ecs table with the label as a pseudo-pathway
  pw_ecs <- data.table(pathway = ec_label, ec = direct_ecs)
  matched_pws <- ec_label

} else {
  # Pathway mode: resolve pathway -> reactions -> ECs via mapping files
  cat("\n[1/5] Resolving pathways to EC numbers...\n")

  # ── Locate PICRUSt2 mapping files ──
  find_picrust2_mapdir <- function() {
    for (py in c("python3", "python")) {
      cmd <- sprintf('%s -c "import picrust2, os; print(os.path.join(os.path.dirname(picrust2.__file__), \'default_files\', \'pathway_mapfiles\'))" 2>/dev/null', py)
      result <- suppressWarnings(system(cmd, intern = TRUE))
      if (length(result) == 1L && dir.exists(result)) return(result)
    }
    conda_base <- Sys.getenv("CONDA_PREFIX", unset = NA)
    if (is.na(conda_base)) conda_base <- Sys.getenv("HOME")
    search_dirs <- c(
      file.path(conda_base, "lib"),
      file.path(Sys.getenv("HOME"), "miniconda3/envs"),
      file.path(Sys.getenv("HOME"), "anaconda3/envs"),
      "/opt/conda/envs"
    )
    for (d in search_dirs) {
      hits <- list.files(d, pattern = "pathway_mapfiles$",
                         recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
      hits <- hits[dir.exists(hits)]
      if (length(hits) > 0L) return(hits[1L])
    }
    NULL
  }

  resolve_mapping_files <- function(args) {
    path2rxn <- args[["path2rxn"]]
    rxn2ec   <- args[["rxn2ec"]]
    if (is.null(path2rxn) || is.null(rxn2ec)) {
      mapdir <- find_picrust2_mapdir()
      if (is.null(mapdir)) {
        stop("Cannot auto-detect PICRUSt2 mapping files. ",
             "Provide --path2rxn and --rxn2ec manually.", call. = FALSE)
      }
      cat("Auto-detected PICRUSt2 mapping files at:\n  ", mapdir, "\n")
      if (is.null(path2rxn))
        path2rxn <- file.path(mapdir, "metacyc_path2rxn_struc_filt_pro.txt")
      if (is.null(rxn2ec))
        rxn2ec <- file.path(mapdir, "metacyc_rxn_to_level4ec.tsv")
    }
    for (f in c(path2rxn, rxn2ec)) {
      if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
    }
    list(path2rxn = path2rxn, rxn2ec = rxn2ec)
  }

  maps <- resolve_mapping_files(args)

  # Parse pathway-to-reaction file (structured MinPath format)
  path2rxn_raw <- readLines(maps$path2rxn)
  path2rxn_dt  <- rbindlist(lapply(path2rxn_raw, function(line) {
    parts <- strsplit(line, "\t")[[1L]]
    if (length(parts) < 2L) return(NULL)
    pw <- trimws(parts[1L])
    rxn_str <- parts[2L]
    tokens <- unlist(regmatches(rxn_str, gregexpr("[A-Za-z0-9_.+-]+[-][A-Za-z0-9_.+-]*|[A-Za-z0-9_.]+", rxn_str)))
    tokens <- tokens[!tokens %in% c("and", "or")]
    tokens <- sub("^-", "", tokens)
    if (length(tokens) == 0L) return(NULL)
    data.table(pathway = pw, reaction = tokens)
  }))

  # Load reaction-to-EC mapping (variable columns: one reaction can map to multiple ECs)
  rxn2ec_lines <- readLines(maps$rxn2ec)
  rxn2ec <- rbindlist(lapply(rxn2ec_lines, function(l) {
    parts <- strsplit(l, "\t")[[1L]]
    if (length(parts) < 2L) return(NULL)
    data.table(reaction = parts[1L], ec = sub("^EC:", "", parts[-1L]))
  }))

  # Filter to requested pathways
  pw_rxns <- path2rxn_dt[pathway %in% pathway_ids]
  if (nrow(pw_rxns) == 0L) {
    found_pws <- unique(path2rxn_dt$pathway)
    cat("No matching pathways found.\n")
    cat("Available pathways (first 20):\n")
    cat(paste(" ", head(found_pws, 20)), sep = "\n")
    stop("None of the requested pathways matched: ",
         paste(pathway_ids, collapse = ", "), call. = FALSE)
  }

  matched_pws <- intersect(pathway_ids, unique(pw_rxns$pathway))
  missed_pws  <- setdiff(pathway_ids, matched_pws)
  if (length(missed_pws) > 0L) {
    cat("  WARNING: pathways not found in mapping file: ",
        paste(missed_pws, collapse = ", "), "\n")
  }

  # Join to get ECs
  pw_ecs <- merge(pw_rxns, rxn2ec, by = "reaction", all.x = FALSE)
  pw_ecs <- pw_ecs[!is.na(ec) & ec != ""]

  for (pw in matched_pws) {
    n_rxn <- uniqueN(pw_rxns[pathway == pw, reaction])
    n_ec  <- uniqueN(pw_ecs[pathway == pw, ec])
    cat(sprintf("  %s: %d reactions -> %d EC numbers\n", pw, n_rxn, n_ec))
  }
}

# ── Step 2: Parse taxonomy ───────────────────────────────────────────────────

cat("\n[2/5] Parsing taxonomy...\n")

tax_raw <- fread(tax_path, sep = "\t")

# Standardize column names
id_col <- names(tax_raw)[1L]
setnames(tax_raw, id_col, "taxon_id")

if ("Taxon" %in% names(tax_raw)) {
  taxon_col <- "Taxon"
} else if ("taxon" %in% names(tax_raw)) {
  taxon_col <- "taxon"
} else {
  taxon_col <- names(tax_raw)[2L]
  cat("  Using column '", taxon_col, "' as taxonomy string\n")
}

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

parse_qiime2_taxonomy <- function(taxon_strings) {
  # Split on semicolon
  split_tax <- strsplit(taxon_strings, ";\\s*")
  # Strip any rank prefix (d__/k__/p__/D_0__/D_1__/etc.)
  parsed <- lapply(split_tax, function(x) {
    cleaned <- sub("^[dkpcofgsDKPCOFGS](_[0-9]*)?__", "", trimws(x))
    # Pad to 7 levels
    length(cleaned) <- 7L
    cleaned
  })
  dt <- as.data.table(do.call(rbind, parsed))
  setnames(dt, rank_names[seq_len(ncol(dt))])
  dt
}

tax_parsed <- parse_qiime2_taxonomy(tax_raw[[taxon_col]])
tax <- cbind(tax_raw[, .(taxon_id)], tax_parsed)

# Replace empty/NA with "Unclassified"
for (col in rank_names) {
  if (col %in% names(tax)) {
    tax[is.na(get(col)) | get(col) == "" | get(col) == "__", (col) := "Unclassified"]
  }
}

cat(sprintf("  %d taxa loaded (%d classified to Family, %d to Genus)\n",
            nrow(tax),
            sum(tax$Family != "Unclassified"),
            sum(tax$Genus != "Unclassified")))

# ── Step 3: Load EC contributor table ────────────────────────────────────────

cat("\n[3/5] Loading EC contributor table (this may take a minute)...\n")

if (!file.exists(ec_contrib_path)) {
  stop("EC contributor file not found: ", ec_contrib_path, call. = FALSE)
}

ec_contrib <- fread(ec_contrib_path)
setnames(ec_contrib, "function", "ec")
ec_contrib[, ec := sub("^(ec|ko):", "", ec)]

cat(sprintf("  %s rows, %d samples, %d ECs\n",
            format(nrow(ec_contrib), big.mark = ","),
            uniqueN(ec_contrib$sample),
            uniqueN(ec_contrib$ec)))

# Filter to target ECs
all_target_ecs <- unique(pw_ecs$ec)
ec_sub <- ec_contrib[ec %in% all_target_ecs]
cat(sprintf("  %s rows match target ECs\n", format(nrow(ec_sub), big.mark = ",")))

if (nrow(ec_sub) == 0L) {
  stop("No matching ECs found in contributor file. ",
       "Check that --ec-contrib points to the correct file.", call. = FALSE)
}

# Filter to requested samples if specified
if (!is.null(sample_ids)) {
  ec_sub <- ec_sub[sample %in% sample_ids]
  cat(sprintf("  Filtered to %d samples: %s rows\n",
              length(sample_ids), format(nrow(ec_sub), big.mark = ",")))
}

# ── Step 4: Merge with taxonomy and aggregate ────────────────────────────────

cat("\n[4/5] Merging with taxonomy and aggregating...\n")

# Rename taxon column in contributor to match taxonomy
setnames(ec_sub, "taxon", "taxon_id")

# Join taxonomy
ec_tax <- merge(ec_sub, tax, by = "taxon_id", all.x = TRUE)
# Fill unmatched taxa
for (col in rank_names) {
  if (col %in% names(ec_tax)) {
    ec_tax[is.na(get(col)), (col) := "Unclassified"]
  }
}

# Add pathway label to each row (an EC can belong to multiple pathways)
ec_pw <- merge(ec_tax, pw_ecs[, .(pathway, ec)], by = "ec", allow.cartesian = TRUE)

# Load metadata and assign groups if requested
has_groups <- FALSE
if (!is.null(meta_path) && !is.null(group_col)) {
  meta <- fread(meta_path)
  meta_id_col <- names(meta)[1L]
  setnames(meta, meta_id_col, "sample")
  if (!group_col %in% names(meta)) {
    cat("  WARNING: column '", group_col, "' not in metadata. Skipping group analysis.\n")
  } else {
    ec_pw <- merge(ec_pw, meta[, .(sample, group = get(group_col))], by = "sample", all.x = TRUE)
    ec_pw[is.na(group), group := "Unknown"]
    has_groups <- TRUE
    cat(sprintf("  Groups defined: %s\n", paste(sort(unique(ec_pw$group)), collapse = ", ")))
  }
}

# ── Create output directory ──────────────────────────────────────────────────

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Step 5: Write output tables ──────────────────────────────────────────────

cat("\n[5/5] Writing output tables...\n")

aggregate_and_write <- function(dt, pathway_id, out_dir, has_groups, top_n) {

  pw_data <- dt[pathway == pathway_id]
  if (nrow(pw_data) == 0L) return(invisible(NULL))
  safe_id <- gsub("[^A-Za-z0-9_-]", "_", pathway_id)

  # ── ASV-level ──
  if (has_groups) {
    asv_agg <- pw_data[, .(
      contribution = sum(taxon_function_abun, na.rm = TRUE),
      n_ecs        = uniqueN(ec),
      n_samples    = uniqueN(sample)
    ), by = .(group, taxon_id, Kingdom, Phylum, Class, Order, Family, Genus, Species)]
    setorder(asv_agg, group, -contribution)
  } else {
    asv_agg <- pw_data[, .(
      contribution = sum(taxon_function_abun, na.rm = TRUE),
      n_ecs        = uniqueN(ec),
      n_samples    = uniqueN(sample)
    ), by = .(taxon_id, Kingdom, Phylum, Class, Order, Family, Genus, Species)]
    setorder(asv_agg, -contribution)
  }
  fwrite(asv_agg, file.path(out_dir, paste0(safe_id, "_ASV.tsv")), sep = "\t")

  # ── Genus-level ──
  group_by_genus <- c(if (has_groups) "group", "Phylum", "Class", "Order", "Family", "Genus")
  genus_agg <- pw_data[, .(
    contribution = sum(taxon_function_abun, na.rm = TRUE),
    n_asvs       = uniqueN(taxon_id),
    n_ecs        = uniqueN(ec),
    n_samples    = uniqueN(sample)
  ), by = group_by_genus]
  if (has_groups) {
    genus_agg[, pct := round(100 * contribution / sum(contribution), 2), by = group]
    setorder(genus_agg, group, -contribution)
  } else {
    genus_agg[, pct := round(100 * contribution / sum(contribution), 2)]
    setorder(genus_agg, -contribution)
  }
  fwrite(genus_agg, file.path(out_dir, paste0(safe_id, "_Genus.tsv")), sep = "\t")

  # ── Family-level ──
  group_by_fam <- c(if (has_groups) "group", "Phylum", "Class", "Order", "Family")
  fam_agg <- pw_data[, .(
    contribution = sum(taxon_function_abun, na.rm = TRUE),
    n_asvs       = uniqueN(taxon_id),
    n_ecs        = uniqueN(ec),
    n_samples    = uniqueN(sample),
    top_genus    = {
      g <- .SD[, .(s = sum(taxon_function_abun)), by = Genus]
      g[which.max(s), Genus]
    }
  ), by = group_by_fam]
  if (has_groups) {
    fam_agg[, pct := round(100 * contribution / sum(contribution), 2), by = group]
    setorder(fam_agg, group, -contribution)
  } else {
    fam_agg[, pct := round(100 * contribution / sum(contribution), 2)]
    setorder(fam_agg, -contribution)
  }
  fwrite(fam_agg, file.path(out_dir, paste0(safe_id, "_Family.tsv")), sep = "\t")

  # ── Group enrichment (if groups present, ≥2 groups) ──
  if (has_groups && uniqueN(pw_data$group) >= 2L) {
    fam_wide <- dcast(fam_agg, Phylum + Class + Order + Family ~ group,
                      value.var = "pct", fill = 0, fun.aggregate = sum)
    grp_names <- setdiff(names(fam_wide), c("Phylum", "Class", "Order", "Family"))
    fam_wide[, mean_other := rowMeans(.SD), .SDcols = grp_names]
    fwrite(fam_wide, file.path(out_dir, paste0(safe_id, "_Family_crosstab.tsv")), sep = "\t")
  }

  # ── EC breakdown per pathway ──
  ec_summary <- pw_data[, .(
    total_contribution = sum(taxon_function_abun, na.rm = TRUE),
    n_asvs             = uniqueN(taxon_id),
    top_family         = {
      f <- .SD[, .(s = sum(taxon_function_abun)), by = Family]
      f[which.max(s), Family]
    },
    top_family_pct = {
      f <- .SD[, .(s = sum(taxon_function_abun)), by = Family]
      round(100 * max(f$s) / sum(f$s), 1)
    }
  ), by = .(ec)]
  setorder(ec_summary, -total_contribution)
  fwrite(ec_summary, file.path(out_dir, paste0(safe_id, "_EC_summary.tsv")), sep = "\t")

  # ── Console summary ──
  cat(sprintf("\n  %s:\n", pathway_id))
  cat(sprintf("    %s total contribution across %d ASVs, %d ECs\n",
              format(round(sum(pw_data$taxon_function_abun)), big.mark = ","),
              uniqueN(pw_data$taxon_id), uniqueN(pw_data$ec)))

  # Console: show top classified families (skip "Unclassified" for readability)
  if (has_groups) {
    for (g in sort(unique(fam_agg$group))) {
      cat(sprintf("    [%s] Top families: ", g))
      classified <- fam_agg[group == g & Family != "Unclassified"]
      top <- head(classified, min(5L, top_n))
      unclass_pct <- fam_agg[group == g & Family == "Unclassified", sum(pct)]
      cat(paste(sprintf("%s (%.1f%%)", top$Family, top$pct), collapse = ", "))
      if (length(unclass_pct) > 0L && unclass_pct > 0)
        cat(sprintf(" [+%.1f%% unclassified]", unclass_pct))
      cat("\n")
    }
  } else {
    cat("    Top families: ")
    classified <- fam_agg[Family != "Unclassified"]
    top <- head(classified, min(5L, top_n))
    unclass_pct <- fam_agg[Family == "Unclassified", sum(pct)]
    cat(paste(sprintf("%s (%.1f%%)", top$Family, top$pct), collapse = ", "))
    if (length(unclass_pct) > 0L && unclass_pct > 0)
      cat(sprintf(" [+%.1f%% unclassified]", unclass_pct))
    cat("\n")
  }

  invisible(NULL)
}

for (pw in matched_pws) {
  aggregate_and_write(ec_pw, pw, out_dir, has_groups, top_n)
}

# ── Final summary ───────────────────────────────────────────────────────────

n_files <- length(list.files(out_dir, pattern = "\\.tsv$"))
cat(sprintf("\nDone. %d files written to %s/\n", n_files, out_dir))
cat("Output files per pathway:\n")
cat("  *_ASV.tsv           — Per-ASV contributions with full taxonomy\n")
cat("  *_Genus.tsv         — Aggregated by genus\n")
cat("  *_Family.tsv        — Aggregated by family\n")
cat("  *_EC_summary.tsv    — Per-EC breakdown with top contributing family\n")
if (has_groups) {
  cat("  *_Family_crosstab.tsv — Family % contribution per group (wide format)\n")
}
