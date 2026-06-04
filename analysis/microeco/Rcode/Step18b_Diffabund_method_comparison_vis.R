##
## Cross-method comparison of differential abundance results
## Compares metagenomeSeq, ANCOM-BC2, and DESeq2 (ASV level) with
## LEfSe (Family level) for both site and season contrasts.
##
## Visualisations produced:
##   1. Lollipop/dot plot of top families — site comparison
##   2. Lollipop/dot plot of top families — season comparison
##   3. Side-by-side LEfSe LDA vs count-based method family heatmap — site
##   4. Side-by-side LEfSe LDA vs count-based method family heatmap — season
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
theme_set(theme_bw())

######################################################
output_dir <- "Output/1.Amplicon/Stage6_Diff_abund"
output_dir_comp <- file.path(output_dir, "3.method_comparison")
dir.create(output_dir_comp, recursive = TRUE, showWarnings = FALSE)


######################################################
# helper: extract Family name from full taxonomic path
# e.g. "k__...|f__Nostocaceae|g__...|s__...|ASV_xxx" -> "Nostocaceae"
extract_family <- function(taxa_string){
  fam <- gsub(".*f__([^|]+)\\|.*", "\\1", taxa_string)
  fam[fam == taxa_string] <- NA          # no f__ match
  fam <- gsub("^$", NA, fam)
  fam
}

# helper: clean family label for display
clean_family <- function(x) gsub("^f__", "", x)


######################################################
# load and parse count-based method results

load_metagenomeseq <- function(path, sig_levels = c("*","**","***")){
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$Significance %in% sig_levels, ]
  df$Family    <- extract_family(df$Taxa)
  df$direction <- ifelse(df$logFC > 0, "Group1", "Group2")
  df[, c("Taxa","Family","logFC","P.adj","Significance","Comparison","direction")]
}

load_deseq2 <- function(path, sig_levels = c("*","**","***")){
  df <- read.csv(path, stringsAsFactors = FALSE)
  df <- df[df$Significance %in% sig_levels, ]
  df$Family    <- extract_family(df$Taxa)
  df$direction <- ifelse(df$log2FoldChange > 0, "Group1", "Group2")
  df[, c("Taxa","Family","log2FoldChange","P.adj","Significance","Comparison","direction")]
}

load_ancombc2 <- function(path){
  df <- read.csv(path, stringsAsFactors = FALSE)
  # keep only the siteEsteroZacatecas / seasonSeptember factor rows
  # (remove Intercept rows which duplicate every ASV)
  df <- df[!grepl("Intercept", df$Factors), ]
  df <- df[df$diff == TRUE, ]
  df$Family    <- extract_family(df$Taxa)
  df$direction <- ifelse(df$lfc > 0, "Group1", "Group2")
  df[, c("Taxa","Family","lfc","P.adj","Significance","direction")]
}


######################################################
# summarise results at Family level
# returns: family, n_sig (number of significant ASVs), dominant direction, mean effect

family_summary <- function(df, effect_col, method_name){
  df <- df[!is.na(df$Family), ]
  df %>%
    group_by(Family) %>%
    summarise(
      n_sig     = n(),
      mean_eff  = mean(abs(.data[[effect_col]]), na.rm = TRUE),
      direction = names(which.max(table(direction))),
      .groups = "drop"
    ) %>%
    mutate(method = method_name)
}


######################################################
##  SITE COMPARISON
######################################################

mgs_site  <- load_metagenomeseq(file.path(output_dir, "Diff_abund_ASV_site_metagenomeSeq.csv"))
ds2_site  <- load_deseq2(file.path(output_dir, "Diff_abund_ASV_site_DESeq2.csv"))
anc_site  <- load_ancombc2(file.path(output_dir, "Diff_abund_ASV_site_ancombc2.csv"))
lef_site  <- read.csv(file.path(output_dir, "Diff_abund_Family_site_lefse.csv"),
                      stringsAsFactors = FALSE)

# parse direction from LEfSe Group column
lef_site$Family    <- clean_family(gsub(".*f__([^|]+).*", "\\1", lef_site$Taxa))
lef_site$direction <- lef_site$Group
lef_site$LDA_signed <- ifelse(lef_site$Group == "EsteroZacatecas",
                               lef_site$LDA, -lef_site$LDA)

# summarise count-based methods at family level
fam_mgs <- family_summary(mgs_site, "logFC",          "metagenomeSeq")
fam_ds2 <- family_summary(ds2_site, "log2FoldChange", "DESeq2")
fam_anc <- family_summary(anc_site, "lfc",            "ANCOM-BC2")

# use direction: positive = EsteroZacatecas enriched (comparison is Cola - Estero for mgs)
# metagenomeSeq: Comparison "ColaDeBallena - EsteroZacatecas"
#   logFC > 0 means higher in ColaDeBallena (Group1 = Cola)
# DESeq2: Comparison "EsteroZacatecas - ColaDeBallena"
#   log2FC > 0 means higher in EsteroZacatecas (Group1 = Estero)
# Harmonise so direction = site enriched
fam_mgs$site_enriched <- ifelse(fam_mgs$direction == "Group1",
                                "EsteroZacatecas", "ColaDeBallena")
fam_ds2$site_enriched <- ifelse(fam_ds2$direction == "Group1",
                                 "EsteroZacatecas", "ColaDeBallena")
fam_anc$site_enriched <- ifelse(fam_anc$direction == "Group1",
                                 "EsteroZacatecas", "ColaDeBallena")

# combine all methods
all_fam_site <- bind_rows(
  fam_mgs[, c("Family","n_sig","mean_eff","site_enriched","method")],
  fam_ds2[, c("Family","n_sig","mean_eff","site_enriched","method")],
  fam_anc[, c("Family","n_sig","mean_eff","site_enriched","method")]
)

# select families present in at least 2 methods or with n_sig >= 3
top_families_site <- all_fam_site %>%
  group_by(Family) %>%
  summarise(total_sig = sum(n_sig), n_methods = n_distinct(method), .groups = "drop") %>%
  filter(n_methods >= 2 | total_sig >= 3) %>%
  arrange(desc(total_sig)) %>%
  slice_head(n = 25) %>%
  pull(Family)

# also ensure LEfSe significant families are included
lefse_sig_families <- lef_site$Family[lef_site$Significance %in% c("*","**","***")]
top_families_site  <- unique(c(lefse_sig_families, top_families_site))

plot_dat_site <- all_fam_site %>%
  filter(Family %in% top_families_site) %>%
  mutate(
    n_sig_display = pmin(n_sig, 30),  # cap bubble size
    Family = factor(Family, levels = rev(sort(unique(Family))))
  )


######################################################
# Figure 1a: bubble plot — site, count-based methods
pal_site <- c("EsteroZacatecas" = "#1b7837", "ColaDeBallena" = "#762a83")

p_site_count <- ggplot(plot_dat_site,
    aes(x = method, y = Family, size = n_sig_display,
        colour = site_enriched, fill = site_enriched)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(name = "Significant\nASVs", range = c(2, 12),
                        breaks = c(1, 5, 15, 30)) +
  scale_colour_manual(values = pal_site, name = "Enriched in") +
  scale_fill_manual(values   = pal_site, name = "Enriched in") +
  labs(x = NULL, y = NULL,
       title = "Site (EsteroZacatecas vs ColaDeBallena)",
       subtitle = "Significant ASVs per family (count-based methods)") +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 11),
        plot.title  = element_text(size = 12, face = "bold"),
        legend.position = "right")

# Figure 1b: LEfSe LDA bar for site
lef_site_sig <- lef_site[lef_site$Significance %in% c("*","**","***"), ]
lef_site_sig$Family <- factor(lef_site_sig$Family,
                               levels = lef_site_sig$Family[order(lef_site_sig$LDA_signed)])

p_site_lefse <- ggplot(lef_site_sig,
    aes(x = LDA_signed, y = Family, fill = direction)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = pal_site, name = "Enriched in") +
  labs(x = "LDA score (signed)", y = NULL,
       title = "LEfSe — Family level (site)",
       subtitle = "Kruskal-Wallis + LDA score ≥ 2") +
  theme(axis.text.y  = element_text(size = 9),
        plot.title   = element_text(size = 12, face = "bold"),
        legend.position = "none")

g_site <- ggarrange(p_site_lefse, p_site_count, ncol = 2,
                    widths = c(1, 1.6), common.legend = TRUE, legend = "right")
png(file.path(output_dir_comp, "compare_site_LEfSe_vs_countmethods.png"),
    width = 3600, height = 2000, res = 300)
print(g_site)
dev.off()


######################################################
##  SEASON COMPARISON
######################################################

mgs_sea  <- load_metagenomeseq(file.path(output_dir, "Diff_abund_ASV_season_metagenomeSeq.csv"))
ds2_sea  <- load_deseq2(file.path(output_dir, "Diff_abund_ASV_season_DESeq2.csv"))
lef_sea  <- read.csv(file.path(output_dir, "Diff_abund_Family_season_lefse.csv"),
                     stringsAsFactors = FALSE)

# ANCOM-BC2 for season (Factors column = "seasonSeptember")
anc_sea  <- load_ancombc2(file.path(output_dir, "Diff_abund_ASV_season_ancombc2.csv"))

lef_sea$Family    <- clean_family(gsub(".*f__([^|]+).*", "\\1", lef_sea$Taxa))
lef_sea$direction <- lef_sea$Group
lef_sea$LDA_signed <- ifelse(lef_sea$Group == "September",
                              lef_sea$LDA, -lef_sea$LDA)

# Harmonise directions for season
# metagenomeSeq: "February - September" → logFC > 0 means higher in February
# DESeq2: "September - February"       → log2FC > 0 means higher in September
fam_mgs_sea <- family_summary(mgs_sea, "logFC",          "metagenomeSeq")
fam_ds2_sea <- family_summary(ds2_sea, "log2FoldChange", "DESeq2")
fam_anc_sea <- family_summary(anc_sea, "lfc",            "ANCOM-BC2")

fam_mgs_sea$season_enriched <- ifelse(fam_mgs_sea$direction == "Group1",
                                      "September", "February")
fam_ds2_sea$season_enriched <- ifelse(fam_ds2_sea$direction == "Group1",
                                       "September", "February")
fam_anc_sea$season_enriched <- ifelse(fam_anc_sea$direction == "Group1",
                                       "September", "February")

all_fam_sea <- bind_rows(
  fam_mgs_sea[, c("Family","n_sig","mean_eff","season_enriched","method")],
  fam_ds2_sea[, c("Family","n_sig","mean_eff","season_enriched","method")],
  fam_anc_sea[, c("Family","n_sig","mean_eff","season_enriched","method")]
)

top_families_sea <- all_fam_sea %>%
  group_by(Family) %>%
  summarise(total_sig = sum(n_sig), n_methods = n_distinct(method), .groups = "drop") %>%
  filter(n_methods >= 2 | total_sig >= 3) %>%
  arrange(desc(total_sig)) %>%
  slice_head(n = 25) %>%
  pull(Family)

lefse_sig_sea    <- lef_sea$Family[lef_sea$Significance %in% c("*","**","***")]
top_families_sea <- unique(c(lefse_sig_sea, top_families_sea))

plot_dat_sea <- all_fam_sea %>%
  filter(Family %in% top_families_sea) %>%
  mutate(
    n_sig_display = pmin(n_sig, 30),
    Family = factor(Family, levels = rev(sort(unique(Family))))
  )

pal_sea <- c("September" = "#d95f02", "February" = "#1f78b4")

p_sea_count <- ggplot(plot_dat_sea,
    aes(x = method, y = Family, size = n_sig_display,
        colour = season_enriched, fill = season_enriched)) +
  geom_point(shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(name = "Significant\nASVs", range = c(2, 12),
                        breaks = c(1, 5, 15, 30)) +
  scale_colour_manual(values = pal_sea, name = "Enriched in") +
  scale_fill_manual(values   = pal_sea, name = "Enriched in") +
  labs(x = NULL, y = NULL,
       title = "Season (February vs September)",
       subtitle = "Significant ASVs per family (count-based methods)") +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 11),
        plot.title  = element_text(size = 12, face = "bold"),
        legend.position = "right")

lef_sea_sig <- lef_sea[lef_sea$Significance %in% c("*","**","***"), ]
lef_sea_sig$Family <- factor(lef_sea_sig$Family,
                              levels = lef_sea_sig$Family[order(lef_sea_sig$LDA_signed)])

p_sea_lefse <- ggplot(lef_sea_sig,
    aes(x = LDA_signed, y = Family, fill = direction)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = pal_sea, name = "Enriched in") +
  labs(x = "LDA score (signed)", y = NULL,
       title = "LEfSe — Family level (season)",
       subtitle = "Kruskal-Wallis + LDA score ≥ 2") +
  theme(axis.text.y  = element_text(size = 9),
        plot.title   = element_text(size = 12, face = "bold"),
        legend.position = "none")

g_sea <- ggarrange(p_sea_lefse, p_sea_count, ncol = 2,
                   widths = c(1, 1.6), common.legend = TRUE, legend = "right")
png(file.path(output_dir_comp, "compare_season_LEfSe_vs_countmethods.png"),
    width = 3600, height = 2000, res = 300)
print(g_sea)
dev.off()


######################################################
##  COMBINED FIGURE — top site ASVs across methods
##  Focus: show the Nostocaceae / Rivularia signal
######################################################

# select ASVs significant in both metagenomeSeq AND DESeq2 (consensus)
sig_both_site <- intersect(mgs_site$Taxa, ds2_site$Taxa)

cat("ASVs significant in both metagenomeSeq and DESeq2 (site):", length(sig_both_site), "\n")

if(length(sig_both_site) > 0){
  # build long table with effect sizes
  mgs_consensus <- mgs_site[mgs_site$Taxa %in% sig_both_site,
                             c("Taxa","Family","logFC","Significance")]
  mgs_consensus$method <- "metagenomeSeq"
  # logFC > 0 = enriched in EsteroZacatecas (second group in "Cola - Estero" comparison)
  mgs_consensus$effect <- mgs_consensus$logFC

  ds2_consensus <- ds2_site[ds2_site$Taxa %in% sig_both_site,
                             c("Taxa","Family","log2FoldChange","Significance")]
  ds2_consensus$method <- "DESeq2"
  ds2_consensus$effect <- ds2_consensus$log2FoldChange

  anc_in_both  <- anc_site[anc_site$Taxa %in% sig_both_site, c("Taxa","Family","lfc","Significance")]
  anc_in_both$method <- "ANCOM-BC2"
  anc_in_both$effect <- anc_in_both$lfc

  long_consensus <- bind_rows(
    mgs_consensus[, c("Taxa","Family","effect","Significance","method")],
    ds2_consensus[, c("Taxa","Family","effect","Significance","method")],
    anc_in_both[,  c("Taxa","Family","effect","Significance","method")]
  )

  # extract short ASV label
  long_consensus$ASV_id <- gsub(".*ASV_(\\d+)$", "ASV_\\1", long_consensus$Taxa)
  long_consensus$label  <- paste0(long_consensus$ASV_id, " [",
                                   long_consensus$Family, "]")

  # keep top 30 by maximum absolute effect across methods
  top_asvs <- long_consensus %>%
    group_by(Taxa) %>%
    summarise(max_eff = max(abs(effect), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_eff)) %>%
    slice_head(n = 30) %>%
    pull(Taxa)

  plot_consensus <- long_consensus %>%
    filter(Taxa %in% top_asvs) %>%
    mutate(label = factor(label, levels = unique(label[order(effect[method == "DESeq2"])])))

  p_consensus <- ggplot(plot_consensus,
      aes(x = effect, y = label, colour = method, shape = method)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = c("metagenomeSeq" = "#1b9e77",
                                   "DESeq2"        = "#d95f02",
                                   "ANCOM-BC2"     = "#7570b3"),
                        name = "Method") +
    scale_shape_manual(values = c("metagenomeSeq" = 16,
                                  "DESeq2"        = 17,
                                  "ANCOM-BC2"     = 15),
                       name = "Method") +
    labs(x = "Effect size (positive = enriched in EsteroZacatecas)",
         y = NULL,
         title = "Consensus significant ASVs — site comparison",
         subtitle = paste0("ASVs significant in both metagenomeSeq and DESeq2 (n=",
                           length(sig_both_site), "); top 30 by effect size")) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom")

  png(file.path(output_dir_comp, "consensus_ASVs_site_effectsize.png"),
      width = 3000, height = 3200, res = 300)
  print(p_consensus)
  dev.off()
}


######################################################
##  COMBINED FIGURE — top season ASVs across methods

sig_both_sea <- intersect(mgs_sea$Taxa, ds2_sea$Taxa)

cat("ASVs significant in both metagenomeSeq and DESeq2 (season):", length(sig_both_sea), "\n")

if(length(sig_both_sea) > 0){
  mgs_sea_con <- mgs_sea[mgs_sea$Taxa %in% sig_both_sea,
                          c("Taxa","Family","logFC","Significance")]
  mgs_sea_con$method <- "metagenomeSeq"
  # logFC > 0 = enriched in September (second group in "February - September" comparison)
  mgs_sea_con$effect <- mgs_sea_con$logFC

  ds2_sea_con <- ds2_sea[ds2_sea$Taxa %in% sig_both_sea,
                          c("Taxa","Family","log2FoldChange","Significance")]
  ds2_sea_con$method <- "DESeq2"
  # DESeq2 comparison "September - February": log2FC > 0 = enriched in September
  ds2_sea_con$effect <- ds2_sea_con$log2FoldChange

  anc_sea_con <- anc_sea[anc_sea$Taxa %in% sig_both_sea,
                          c("Taxa","Family","lfc","Significance")]
  anc_sea_con$method <- "ANCOM-BC2"
  # lfc for seasonSeptember factor: > 0 = enriched in September
  anc_sea_con$effect <- anc_sea_con$lfc

  long_sea_con <- bind_rows(
    mgs_sea_con[, c("Taxa","Family","effect","Significance","method")],
    ds2_sea_con[, c("Taxa","Family","effect","Significance","method")],
    anc_sea_con[, c("Taxa","Family","effect","Significance","method")]
  )

  long_sea_con$ASV_id <- gsub(".*ASV_(\\d+)$", "ASV_\\1", long_sea_con$Taxa)
  long_sea_con$label  <- paste0(long_sea_con$ASV_id, " [", long_sea_con$Family, "]")

  top_asvs_sea <- long_sea_con %>%
    group_by(Taxa) %>%
    summarise(max_eff = max(abs(effect), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_eff)) %>%
    slice_head(n = 30) %>%
    pull(Taxa)

  plot_sea_con <- long_sea_con %>%
    filter(Taxa %in% top_asvs_sea) %>%
    mutate(label = factor(label,
                          levels = unique(label[order(effect[method == "DESeq2"])])))

  p_sea_consensus <- ggplot(plot_sea_con,
      aes(x = effect, y = label, colour = method, shape = method)) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = c("metagenomeSeq" = "#1b9e77",
                                   "DESeq2"        = "#d95f02",
                                   "ANCOM-BC2"     = "#7570b3"),
                        name = "Method") +
    scale_shape_manual(values = c("metagenomeSeq" = 16,
                                  "DESeq2"        = 17,
                                  "ANCOM-BC2"     = 15),
                       name = "Method") +
    labs(x = "Effect size (positive = enriched in September)",
         y = NULL,
         title = "Consensus significant ASVs — season comparison",
         subtitle = paste0("ASVs significant in both metagenomeSeq and DESeq2 (n=",
                           length(sig_both_sea), "); top 30 by effect size")) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom")

  png(file.path(output_dir_comp, "consensus_ASVs_season_effectsize.png"),
      width = 3000, height = 3200, res = 300)
  print(p_sea_consensus)
  dev.off()
}


######################################################
##  SUMMARY TABLE — family-level method concordance (site)

fam_site_wide <- all_fam_site %>%
  select(Family, n_sig, site_enriched, method) %>%
  pivot_wider(names_from = method, values_from = c(n_sig, site_enriched),
              names_glue = "{method}_{.value}") %>%
  left_join(
    lef_site[lef_site$Significance %in% c("*","**","***"),
             c("Family","LDA","direction","Significance")] %>%
      rename(LEfSe_LDA = LDA, LEfSe_enriched = direction, LEfSe_sig = Significance),
    by = "Family"
  ) %>%
  arrange(desc(rowSums(select(., starts_with("metagenomeSeq_n") |
                                 starts_with("DESeq2_n") |
                                 starts_with("ANCOM-BC2_n")),
                       na.rm = TRUE)))

write.csv(fam_site_wide,
          file.path(output_dir_comp, "family_concordance_site.csv"),
          row.names = FALSE)

fam_sea_wide <- all_fam_sea %>%
  select(Family, n_sig, season_enriched, method) %>%
  pivot_wider(names_from = method, values_from = c(n_sig, season_enriched),
              names_glue = "{method}_{.value}") %>%
  left_join(
    lef_sea[lef_sea$Significance %in% c("*","**","***"),
            c("Family","LDA","direction","Significance")] %>%
      rename(LEfSe_LDA = LDA, LEfSe_enriched = direction, LEfSe_sig = Significance),
    by = "Family"
  ) %>%
  arrange(desc(rowSums(select(., starts_with("metagenomeSeq_n") |
                                 starts_with("DESeq2_n") |
                                 starts_with("ANCOM-BC2_n")),
                       na.rm = TRUE)))

write.csv(fam_sea_wide,
          file.path(output_dir_comp, "family_concordance_season.csv"),
          row.names = FALSE)

######################################################
##  CUSTOM VOLCANO PLOTS — genus/family labels on significant ASVs
##  Replaces microeco's built-in plot_volcano() which shows raw ASV IDs

# helper: extract genus label, fall back to family if genus is empty
asv_label <- function(taxa_string){
  genus  <- gsub(".*g__([^|]+)\\|.*", "\\1", taxa_string)
  family <- gsub(".*f__([^|]+)\\|.*", "\\1", taxa_string)
  asv    <- gsub(".*ASV_(\\d+)$", "ASV_\\1", taxa_string)
  label  <- ifelse(nchar(genus) > 0 & genus != taxa_string & !grepl("^s__$", genus),
                   genus, family)
  # append ASV id to distinguish multiple ASVs from the same genus
  paste0(label, " (", asv, ")")
}

# colour scheme: grey = not significant, gradient by direction for significant
# y-axis uses unadjusted p-values for proper volcano shape;
# point colour and labels are based on adjusted p-value significance.
# The dashed threshold line is drawn at the unadjusted p corresponding to
# the most lenient p_adj that is still significant (i.e. the largest p_unadj
# among significant features), giving an empirical FDR-based cutoff.
volcano_plot <- function(df, effect_col, padj_col, punadj_col, sig_col,
                         group1_label, group2_label,
                         title_text,
                         sig_levels = c("**","***"),
                         max_label  = 20){
  df$effect   <- df[[effect_col]]
  # y-axis: unadjusted p for proper volcano shape
  df$neglogp  <- -log10(pmax(df[[punadj_col]], 1e-300))
  df$sig      <- df[[sig_col]] %in% sig_levels
  df$label    <- asv_label(df$Taxa)

  # direction only for significant ASVs; positive = group2 convention
  df$direction <- ifelse(!df$sig, "ns",
                   ifelse(df$effect > 0, group2_label, group1_label))

  pal <- c("ns" = "grey70")
  pal[group1_label] <- "#762a83"
  pal[group2_label] <- "#1b7837"

  # threshold line: largest unadjusted p among adjusted-significant features
  sig_punadj <- df[[punadj_col]][df$sig]
  threshold  <- if(length(sig_punadj) > 0) max(sig_punadj, na.rm = TRUE) else 0.05

  # label top max_label significant points by smallest adjusted p
  df_lab <- df[df$sig, ]
  if(nrow(df_lab) > max_label){
    df_lab <- df_lab[order(df_lab[[padj_col]]), ][1:max_label, ]
  }

  ggplot(df, aes(x = effect, y = neglogp, colour = direction)) +
    geom_point(size = 1.2, alpha = 0.6) +
    ggrepel::geom_text_repel(data = df_lab,
      aes(label = label), size = 2.8, max.overlaps = 20,
      segment.size = 0.3, segment.colour = "grey50") +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_colour_manual(values = pal, name = "Enriched in",
                        breaks = c(group1_label, group2_label)) +
    labs(x = paste0("log fold change  (positive = ", group2_label, ")"),
         y = expression(-log[10](p[unadj])),
         caption = "Dashed line: unadjusted p threshold corresponding to FDR-significant features",
         title = title_text) +
    theme(legend.position = "right",
          plot.title   = element_text(size = 12, face = "bold"),
          plot.caption = element_text(size = 7, colour = "grey50"))
}

# install ggrepel if needed
if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# --- site: metagenomeSeq ---
# logFC > 0 = enriched in EsteroZacatecas (confirmed convention)
site_mgs_df <- read.csv(file.path(output_dir, "Diff_abund_ASV_site_metagenomeSeq.csv"),
                         stringsAsFactors = FALSE)
p_vol_site_mgs <- volcano_plot(site_mgs_df,
  effect_col = "logFC", padj_col = "P.adj", punadj_col = "P.unadj",
  sig_col = "Significance",
  group1_label = "ColaDeBallena", group2_label = "EsteroZacatecas",
  title_text = "metagenomeSeq — site (ASV level)")
png(file.path(output_dir_comp, "volcano_site_metagenomeSeq.png"),
    width = 2800, height = 2400, res = 300)
print(p_vol_site_mgs)
dev.off()

# --- site: DESeq2 ---
# log2FoldChange > 0 = enriched in EsteroZacatecas ("EsteroZacatecas - ColaDeBallena")
site_ds2_df <- read.csv(file.path(output_dir, "Diff_abund_ASV_site_DESeq2.csv"),
                         stringsAsFactors = FALSE)
p_vol_site_ds2 <- volcano_plot(site_ds2_df,
  effect_col = "log2FoldChange", padj_col = "P.adj", punadj_col = "P.unadj",
  sig_col = "Significance",
  group1_label = "ColaDeBallena", group2_label = "EsteroZacatecas",
  title_text = "DESeq2 — site (ASV level)")
png(file.path(output_dir_comp, "volcano_site_DESeq2.png"),
    width = 2800, height = 2400, res = 300)
print(p_vol_site_ds2)
dev.off()

# --- season: metagenomeSeq ---
# logFC > 0 = enriched in September ("February - September" comparison)
sea_mgs_df <- read.csv(file.path(output_dir, "Diff_abund_ASV_season_metagenomeSeq.csv"),
                        stringsAsFactors = FALSE)
p_vol_sea_mgs <- volcano_plot(sea_mgs_df,
  effect_col = "logFC", padj_col = "P.adj", punadj_col = "P.unadj",
  sig_col = "Significance",
  group1_label = "February", group2_label = "September",
  title_text = "metagenomeSeq — season (ASV level)")
png(file.path(output_dir_comp, "volcano_season_metagenomeSeq.png"),
    width = 2800, height = 2400, res = 300)
print(p_vol_sea_mgs)
dev.off()

# --- season: DESeq2 ---
sea_ds2_df <- read.csv(file.path(output_dir, "Diff_abund_ASV_season_DESeq2.csv"),
                        stringsAsFactors = FALSE)
# confirm direction from Comparison column
cat("DESeq2 season Comparison:", unique(sea_ds2_df$Comparison), "\n")
# If "September - February": log2FC > 0 = September enriched  (group2 = September)
# If "February - September": log2FC > 0 = February enriched  (group2 = February — swap labels below)
p_vol_sea_ds2 <- volcano_plot(sea_ds2_df,
  effect_col = "log2FoldChange", padj_col = "P.adj", punadj_col = "P.unadj",
  sig_col = "Significance",
  group1_label = "February", group2_label = "September",
  title_text = "DESeq2 — season (ASV level)")
png(file.path(output_dir_comp, "volcano_season_DESeq2.png"),
    width = 2800, height = 2400, res = 300)
print(p_vol_sea_ds2)
dev.off()


cat("Done. Outputs saved to:", output_dir_comp, "\n")
cat("Files produced:\n")
cat(" - compare_site_LEfSe_vs_countmethods.png\n")
cat(" - compare_season_LEfSe_vs_countmethods.png\n")
cat(" - consensus_ASVs_site_effectsize.png\n")
cat(" - consensus_ASVs_season_effectsize.png\n")
cat(" - volcano_site_metagenomeSeq.png\n")
cat(" - volcano_site_DESeq2.png\n")
cat(" - volcano_season_metagenomeSeq.png\n")
cat(" - volcano_season_DESeq2.png\n")
cat(" - family_concordance_site.csv\n")
cat(" - family_concordance_season.csv\n")


##### CHECK SEVERAL ASVs for Rivularia_PCC1667
library(microeco)

load(file.path(output_dir, "tmp_microtable_preproc.RData"))
tmp_microtable$cal_abund(rel = TRUE)

# find all Rivularia ASVs
riv_asvs <- rownames(tmp_microtable$tax_table)[
  grepl("Rivularia", tmp_microtable$tax_table$Genus)
]
cat("Rivularia ASVs in dataset:", length(riv_asvs), "\n")
cat(riv_asvs, sep = "\n")

# combined relative abundance per sample
riv_abund <- colSums(tmp_microtable$otu_table[riv_asvs, , drop = FALSE])
riv_df <- data.frame(
  Sample    = names(riv_abund),
  Abundance = as.numeric(riv_abund)
)
riv_df <- merge(riv_df, tmp_microtable$sample_table, by.x = "Sample", by.y
                = "row.names")

# mean by site and season
aggregate(Abundance ~ site, data = riv_df, FUN = mean)
aggregate(Abundance ~ season, data = riv_df, FUN = mean)
aggregate(Abundance ~ site + season, data = riv_df, FUN = mean)

# 1. Confirm it is truly zero at ColaDeBallena (not just near-zero)
riv_per_sample <- data.frame(
  t(tmp_microtable$otu_table[riv_asvs, , drop = FALSE])
)
riv_per_sample$site   <- tmp_microtable$sample_table$site
riv_per_sample$season <- tmp_microtable$sample_table$season
riv_per_sample$total  <- rowSums(riv_per_sample[, riv_asvs])

# any non-zero at ColaDeBallena?
riv_per_sample[riv_per_sample$site == "ColaDeBallena" &
                 riv_per_sample$total > 0, ]

# 2. Relative abundance at EsteroZacatecas to gauge ecological importance
tmp_microtable$cal_abund(rel = TRUE)
# total relative abundance per sample (from abund_list at ASV level)
riv_rel <- colSums(
  tmp_microtable$abund_list$ASV[riv_asvs, , drop = FALSE]
) * 100  # convert to %
riv_rel_df <- data.frame(
  Sample    = names(riv_rel),
  RelAbund  = as.numeric(riv_rel)
)
riv_rel_df <- merge(riv_rel_df, tmp_microtable$sample_table,
                    by.x = "Sample", by.y = "row.names")
aggregate(RelAbund ~ site + season, data = riv_rel_df, FUN = mean)