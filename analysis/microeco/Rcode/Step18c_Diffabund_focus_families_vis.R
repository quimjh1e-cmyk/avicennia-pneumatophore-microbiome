##
## Focused method-comparison plots: consensus + top-LDA families
## Produces two separate PNGs:
##   focus_site_LEfSe_vs_countmethods.png
##   focus_season_LEfSe_vs_countmethods.png
##
## Families shown: consensus across all 3 key methods (metagenomeSeq,
## DESeq2, ANCOM-BC2) + top-5 LDA per group from LEfSe.
## Families ordered by median signed log2FC across 3 count-based methods.
##

library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
theme_set(theme_bw())

output_dir      <- "Output/1.Amplicon/Stage6_Diff_abund"
output_dir_comp <- file.path(output_dir, "3.method_comparison")
dir.create(output_dir_comp, recursive = TRUE, showWarnings = FALSE)


######################################################
# helpers (identical to Step18b)
extract_family <- function(taxa_string){
  fam <- gsub(".*f__([^|]+)\\|.*", "\\1", taxa_string)
  fam[fam == taxa_string] <- NA
  fam <- gsub("^$", NA, fam)
  fam
}
clean_family <- function(x) gsub("^f__", "", x)

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
  df <- df[!grepl("Intercept", df$Factors), ]
  df <- df[df$diff == TRUE, ]
  df$Family    <- extract_family(df$Taxa)
  df$direction <- ifelse(df$lfc > 0, "Group1", "Group2")
  df[, c("Taxa","Family","lfc","P.adj","Significance","direction")]
}

family_summary <- function(df, effect_col, method_name){
  df <- df[!is.na(df$Family), ]
  df %>%
    group_by(Family) %>%
    summarise(
      n_sig      = n(),
      mean_eff   = mean(abs(.data[[effect_col]]), na.rm = TRUE),
      median_lfc = median(.data[[effect_col]], na.rm = TRUE),
      direction  = names(which.max(table(direction))),
      .groups = "drop"
    ) %>%
    mutate(method = method_name)
}


######################################################
## FOCUS FAMILY LISTS
######################################################

# Site: consensus (detected by all 3 key methods) + top LDA per group
focus_site <- c(
  "Paracoccaceae", "Cyclobacteriaceae", "Flavobacteriaceae",
  "Microtrichaceae", "Sphingomonadaceae", "Rhodothermaceae",
  "Saprospiraceae", "Ilumatobacteraceae",
  "Parvularculaceae", "Methyloligellaceae",
  "Nostocaceae", "Xenococcaceae", "Alteromonadaceae", "Neomegalonemataceae",
  "Nocardioidaceae", "Sandaracinaceae", "Cellvibrionaceae",
  "Prolixibacteraceae", "Marinobacteraceae"
)

# Season: consensus (all 3 key methods) + top LDA per group
focus_sea <- c(
  "Parvularculaceae", "Rhodothermaceae", "Sphingomonadaceae",
  "Nitriliruptoraceae", "Euzebyaceae",
  "Paracoccaceae", "Cyclobacteriaceae",
  "Saprospiraceae", "Geodermatophilaceae",
  "Flavobacteriaceae", "Pirellulaceae", "Desulfocapsaceae", "Woeseiaceae",
  "Prolixibacteraceae",
  "Rhizobiaceae", "Nocardioidaceae", "Microtrichaceae"
)


######################################################
##  SITE
######################################################

mgs_site <- load_metagenomeseq(file.path(output_dir, "Diff_abund_ASV_site_metagenomeSeq.csv"))
ds2_site <- load_deseq2(file.path(output_dir, "Diff_abund_ASV_site_DESeq2.csv"))
anc_site <- load_ancombc2(file.path(output_dir, "Diff_abund_ASV_site_ancombc2.csv"))
lef_site <- read.csv(file.path(output_dir, "Diff_abund_Family_site_lefse.csv"),
                     stringsAsFactors = FALSE)

lef_site$Family    <- clean_family(gsub(".*f__([^|]+).*", "\\1", lef_site$Taxa))
lef_site$direction <- lef_site$Group
lef_site$LDA_signed <- ifelse(lef_site$Group == "EsteroZacatecas",
                               lef_site$LDA, -lef_site$LDA)

fam_mgs <- family_summary(mgs_site, "logFC",          "metagenomeSeq")
fam_ds2 <- family_summary(ds2_site, "log2FoldChange", "DESeq2")
fam_anc <- family_summary(anc_site, "lfc",            "ANCOM-BC2")

fam_mgs$site_enriched <- ifelse(fam_mgs$direction == "Group1", "EsteroZacatecas", "ColaDeBallena")
fam_ds2$site_enriched <- ifelse(fam_ds2$direction == "Group1", "EsteroZacatecas", "ColaDeBallena")
fam_anc$site_enriched <- ifelse(fam_anc$direction == "Group1", "EsteroZacatecas", "ColaDeBallena")

all_fam_site <- bind_rows(
  fam_mgs[, c("Family","n_sig","mean_eff","median_lfc","site_enriched","method")],
  fam_ds2[, c("Family","n_sig","mean_eff","median_lfc","site_enriched","method")],
  fam_anc[, c("Family","n_sig","mean_eff","median_lfc","site_enriched","method")]
)

# Bubble plot: order by median signed log2FC across the 3 methods
lef_site_sig <- lef_site %>%
  filter(Significance %in% c("*","**","***"), Family %in% focus_site)

cross_method_order_site <- all_fam_site %>%
  filter(Family %in% focus_site) %>%
  group_by(Family) %>%
  summarise(median_cross = median(median_lfc, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_cross)

missing_site <- setdiff(focus_site, cross_method_order_site$Family)
family_levels_site <- c(missing_site, cross_method_order_site$Family)

# LEfSe panel: keep original LDA-based ordering
lda_order_site <- lef_site_sig %>%
  arrange(LDA_signed) %>%
  pull(Family) %>%
  unique()
lefse_only_site <- setdiff(focus_site, lda_order_site)
family_levels_lefse_site <- c(lefse_only_site, lda_order_site)

lef_site_sig <- lef_site_sig %>%
  mutate(Family = factor(Family, levels = family_levels_lefse_site))

plot_dat_site <- all_fam_site %>%
  filter(Family %in% focus_site) %>%
  mutate(
    n_sig_display = pmin(n_sig, 40),
    Family = factor(Family, levels = family_levels_site)
  )

pal_site <- c("EsteroZacatecas" = "#d95f02", "ColaDeBallena" = "seagreen")

p_site_count <- ggplot(plot_dat_site,
    aes(x = method, y = Family, size = n_sig_display,
        colour = site_enriched, fill = site_enriched)) +
  geom_point(shape = 21, alpha = 0.82, stroke = 0.6) +
  scale_size_continuous(name = "Significant\nASVs", range = c(2, 14),
                        breaks = c(1, 5, 15, 30, 40)) +
  scale_colour_manual(values = pal_site, name = "Enriched in") +
  scale_fill_manual(values   = pal_site, name = "Enriched in") +
  labs(x = NULL, y = NULL,
       title = "Site comparison: EsteroZacatecas vs ColaDeBallena",
       subtitle = "Bubble size = sig. ASVs; ordered by median log2FC") +
  theme(
    axis.text.y   = element_text(size = 9),
    axis.text.x   = element_text(size = 11),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_site_lefse <- ggplot(lef_site_sig,
    aes(x = LDA_signed, y = Family, fill = direction)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = pal_site, name = "Enriched in") +
  labs(x = "LDA score (signed)", y = NULL,
       title = "LEfSe — Family level",
       subtitle = "KW + LDA ≥ 2, ordered by score") +
  theme(
    axis.text.y   = element_text(size = 9),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

g_site <- ggarrange(p_site_lefse, p_site_count,
                    ncol = 2, widths = c(1, 1.7),
                    common.legend = TRUE, legend = "right")

png(file.path(output_dir_comp, "focus_site_LEfSe_vs_countmethods.png"),
    width = 3800, height = 2400, res = 300)
print(g_site)
dev.off()
message("Saved: focus_site_LEfSe_vs_countmethods.png")


######################################################
##  SEASON
######################################################

mgs_sea <- load_metagenomeseq(file.path(output_dir, "Diff_abund_ASV_season_metagenomeSeq.csv"))
ds2_sea <- load_deseq2(file.path(output_dir, "Diff_abund_ASV_season_DESeq2.csv"))
anc_sea <- load_ancombc2(file.path(output_dir, "Diff_abund_ASV_season_ancombc2.csv"))
lef_sea <- read.csv(file.path(output_dir, "Diff_abund_Family_season_lefse.csv"),
                    stringsAsFactors = FALSE)

lef_sea$Family    <- clean_family(gsub(".*f__([^|]+).*", "\\1", lef_sea$Taxa))
lef_sea$direction <- lef_sea$Group
lef_sea$LDA_signed <- ifelse(lef_sea$Group == "September",
                              lef_sea$LDA, -lef_sea$LDA)

fam_mgs_sea <- family_summary(mgs_sea, "logFC",          "metagenomeSeq")
fam_ds2_sea <- family_summary(ds2_sea, "log2FoldChange", "DESeq2")
fam_anc_sea <- family_summary(anc_sea, "lfc",            "ANCOM-BC2")

fam_mgs_sea$season_enriched <- ifelse(fam_mgs_sea$direction == "Group1", "September", "February")
fam_ds2_sea$season_enriched <- ifelse(fam_ds2_sea$direction == "Group1", "September", "February")
fam_anc_sea$season_enriched <- ifelse(fam_anc_sea$direction == "Group1", "September", "February")

all_fam_sea <- bind_rows(
  fam_mgs_sea[, c("Family","n_sig","mean_eff","median_lfc","season_enriched","method")],
  fam_ds2_sea[, c("Family","n_sig","mean_eff","median_lfc","season_enriched","method")],
  fam_anc_sea[, c("Family","n_sig","mean_eff","median_lfc","season_enriched","method")]
)

# Bubble plot: order by median signed log2FC across the 3 methods
lef_sea_sig <- lef_sea %>%
  filter(Significance %in% c("*","**","***"), Family %in% focus_sea)

cross_method_order_sea <- all_fam_sea %>%
  filter(Family %in% focus_sea) %>%
  group_by(Family) %>%
  summarise(median_cross = median(median_lfc, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_cross)

missing_sea <- setdiff(focus_sea, cross_method_order_sea$Family)
family_levels_sea <- c(missing_sea, cross_method_order_sea$Family)

# LEfSe panel: keep original LDA-based ordering
lda_order_sea <- lef_sea_sig %>%
  arrange(LDA_signed) %>%
  pull(Family) %>%
  unique()
lefse_only_sea <- setdiff(focus_sea, lda_order_sea)
family_levels_lefse_sea <- c(lefse_only_sea, lda_order_sea)

lef_sea_sig <- lef_sea_sig %>%
  mutate(Family = factor(Family, levels = family_levels_lefse_sea))

plot_dat_sea <- all_fam_sea %>%
  filter(Family %in% focus_sea) %>%
  mutate(
    n_sig_display = pmin(n_sig, 40),
    Family = factor(Family, levels = family_levels_sea)
  )

pal_sea <- c("September" = "#d95f02", "February" = "seagreen")

p_sea_count <- ggplot(plot_dat_sea,
    aes(x = method, y = Family, size = n_sig_display,
        colour = season_enriched, fill = season_enriched)) +
  geom_point(shape = 21, alpha = 0.82, stroke = 0.6) +
  scale_size_continuous(name = "Significant\nASVs", range = c(2, 14),
                        breaks = c(1, 5, 15, 30, 40)) +
  scale_colour_manual(values = pal_sea, name = "Enriched in") +
  scale_fill_manual(values   = pal_sea, name = "Enriched in") +
  labs(x = NULL, y = NULL,
       title = "Season comparison: February vs September",
       subtitle = "Bubble size = sig. ASVs; ordered by median log2FC") +
  theme(
    axis.text.y   = element_text(size = 9),
    axis.text.x   = element_text(size = 11),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_sea_lefse <- ggplot(lef_sea_sig,
    aes(x = LDA_signed, y = Family, fill = direction)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = pal_sea, name = "Enriched in") +
  labs(x = "LDA score (signed)", y = NULL,
       title = "LEfSe — Family level",
       subtitle = "KW + LDA ≥ 2, ordered by score") +
  theme(
    axis.text.y   = element_text(size = 9),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

g_sea <- ggarrange(p_sea_lefse, p_sea_count,
                   ncol = 2, widths = c(1, 1.7),
                   common.legend = TRUE, legend = "right")

png(file.path(output_dir_comp, "focus_season_LEfSe_vs_countmethods.png"),
    width = 3800, height = 2400, res = 300)
print(g_sea)
dev.off()
message("Saved: focus_season_LEfSe_vs_countmethods.png")