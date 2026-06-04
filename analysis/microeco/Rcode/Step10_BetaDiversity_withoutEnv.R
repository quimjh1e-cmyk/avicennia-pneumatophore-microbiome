##
## Beta diversity analysis — no environmental factors
## Groups: season (February, September) x site (ColaDeBallena, EsteroZacatecas)
##


###########################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage5_BetaDiversity"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
load(input_path)
###########################

# add a combined season x site group column to both microtable objects
amplicon_16S_microtable$sample_table$group <- paste(
	amplicon_16S_microtable$sample_table$season,
	amplicon_16S_microtable$sample_table$site, sep = "_")

amplicon_16S_microtable_rarefy$sample_table$group <- paste(
	amplicon_16S_microtable_rarefy$sample_table$season,
	amplicon_16S_microtable_rarefy$sample_table$site, sep = "_")

# working copies
tmp_raw    <- clone(amplicon_16S_microtable)
tmp_rarefy <- clone(amplicon_16S_microtable_rarefy)

# convert grouping variables to factors so ggplot ellipses work without warnings
tmp_raw$sample_table$season    %<>% factor(levels = c("February", "September"))
tmp_raw$sample_table$site      %<>% factor(levels = c("ColaDeBallena", "EsteroZacatecas"))
tmp_raw$sample_table$group     %<>% factor()
tmp_rarefy$sample_table$season %<>% factor(levels = c("February", "September"))
tmp_rarefy$sample_table$site   %<>% factor(levels = c("ColaDeBallena", "EsteroZacatecas"))
tmp_rarefy$sample_table$group  %<>% factor()


###########################################################################################
# BETA DIVERSITY MATRICES
###########################################################################################

# Bray-Curtis + UniFrac (rarefied; requires rooted phylogenetic tree)
tmp_rarefy$cal_betadiv(method = "bray", unifrac = TRUE)
write.csv(tmp_rarefy$beta_diversity[["bray"]],
	file.path(output_dir, "BetaDiv_rarefy_bray.csv"))
write.csv(tmp_rarefy$beta_diversity$wei_unifrac,
	file.path(output_dir, "BetaDiv_rarefy_weightedUniFrac.csv"))
write.csv(tmp_rarefy$beta_diversity$unwei_unifrac,
	file.path(output_dir, "BetaDiv_rarefy_unweightedUniFrac.csv"))

# Jaccard (binary presence/absence, rarefied)
tmp_rarefy$cal_betadiv(method = "jaccard")
write.csv(tmp_rarefy$beta_diversity[["jaccard"]],
	file.path(output_dir, "BetaDiv_rarefy_jaccard.csv"))

# Aitchison (CLR-based, uses raw non-rarefied data)
tmp_raw$cal_betadiv(method = "aitchison")
write.csv(tmp_raw$beta_diversity[["aitchison"]],
	file.path(output_dir, "BetaDiv_raw_aitchison.csv"))


###########################################################################################
# PCoA
###########################################################################################

# recalculate bray + UniFrac together so both are available for all ordinations below
tmp_rarefy$cal_betadiv(method = "bray", unifrac = TRUE)

# ---- PCoA: Bray-Curtis (rarefied), coloured by season ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "bray")
t1$cal_ordination(method = "PCoA", ncomp = 2)
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_PCoA_Bray_season_Score.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site",
	plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_PCoA_Bray_season.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# ---- PCoA: Bray-Curtis (rarefied), coloured by site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "site", measure = "bray")
t1$cal_ordination(method = "PCoA", ncomp = 2)
g1 <- t1$plot_ordination(plot_color = "site", plot_shape = "season",
	plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_PCoA_Bray_site.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# ---- PCoA: Bray-Curtis (rarefied), coloured by combined group ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "group", measure = "bray")
t1$cal_ordination(method = "PCoA", ncomp = 2)
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_PCoA_Bray_group_Score.csv"))
g1 <- t1$plot_ordination(plot_color = "group", plot_shape = "group",
	plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_PCoA_Bray_group.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# ---- PCoA: weighted UniFrac (rarefied) ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "wei_unifrac")
t1$cal_ordination(method = "PCoA", ncomp = 2)
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_PCoA_wUniFrac_Score.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site",
	plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_PCoA_wUniFrac_season.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# ---- PCoA: Aitchison (raw non-rarefied) ----
t2 <- trans_beta$new(dataset = tmp_raw, group = "season", measure = "aitchison")
t2$cal_ordination(method = "PCoA", ncomp = 2)
write.csv(t2$res_ordination$scores,
	file.path(output_dir, "BetaDiv_raw_PCoA_Aitchison_Score.csv"))
g2 <- t2$plot_ordination(plot_color = "season", plot_shape = "site",
	plot_type = c("point", "ellipse"))
cowplot::save_plot(file.path(output_dir, "BetaDiv_raw_PCoA_Aitchison_season.png"),
	g2, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


###########################################################################################
# NMDS
###########################################################################################

# ---- NMDS: Bray-Curtis (rarefied) ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "bray")
t1$cal_ordination(method = "NMDS")
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_NMDS_Bray_Score.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site",
	plot_type = c("point", "ellipse"),
	NMDS_stress_pos = c(1.1, 1.4), NMDS_stress_text_prefix = "Stress: ")
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_NMDS_Bray_season.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# ---- NMDS: Aitchison (raw) ----
t2 <- trans_beta$new(dataset = tmp_raw, group = "season", measure = "aitchison")
t2$cal_ordination(method = "NMDS")
write.csv(t2$res_ordination$scores,
	file.path(output_dir, "BetaDiv_raw_NMDS_Aitchison_Score.csv"))
g2 <- t2$plot_ordination(plot_color = "season", plot_shape = "site",
	plot_type = c("point", "ellipse"),
	NMDS_stress_pos = c(1.1, 1.3), NMDS_stress_text_prefix = "Stress: ")
cowplot::save_plot(file.path(output_dir, "BetaDiv_raw_NMDS_Aitchison_season.png"),
	g2, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


###########################################################################################
# PCA and DCA at genus level
###########################################################################################

# aggregate to genus level and clean unidentified genera
tmp_rarefy_genus <- tmp_rarefy$merge_taxa(taxa = "Genus")
tmp_rarefy_genus$tax_table %<>% .[.$Genus != "g__", ]
tmp_rarefy_genus$tax_table %<>% .[!grepl("\\d+", .$Genus), ]
tmp_rarefy_genus$tidy_dataset()
rownames(tmp_rarefy_genus$tax_table) <- rownames(tmp_rarefy_genus$otu_table) <-
	gsub("g__", "", tmp_rarefy_genus$tax_table$Genus)

# ---- PCA ----
t1 <- trans_beta$new(dataset = tmp_rarefy_genus)
t1$cal_ordination(method = "PCA", scale_species = TRUE, scale_species_ratio = 1)
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_PCA_Genus_Score.csv"))
write.csv(t1$res_ordination$loading,
	file.path(output_dir, "BetaDiv_rarefy_PCA_Genus_Loading.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site", loading_arrow = TRUE)
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_PCA_Genus.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# ---- DCA ----
t1 <- trans_beta$new(dataset = tmp_rarefy_genus)
t1$cal_ordination(method = "DCA", scale_species = TRUE)
write.csv(t1$res_ordination$scores,
	file.path(output_dir, "BetaDiv_rarefy_DCA_Genus_Score.csv"))
write.csv(t1$res_ordination$loading,
	file.path(output_dir, "BetaDiv_rarefy_DCA_Genus_Loading.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site", loading_arrow = TRUE)
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_DCA_Genus.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)


###########################################################################################
# GROUP DISTANCES
###########################################################################################
            
              ## skipped since there are only two groups ##

# ---- Within-season distances, broken down by site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "bray")
t1$cal_group_distance(within_group = TRUE, by_group = "site")
write.csv(t1$res_group_distance,
	file.path(output_dir, "BetaDiv_rarefy_bray_season_within_bysite.csv"))
t1$cal_group_distance_diff(method = "wilcox")
write.csv(t1$res_group_distance_diff,
	file.path(output_dir, "BetaDiv_rarefy_bray_season_within_bysite_diff.csv"))
g1 <- t1$plot_group_distance(xtext_angle = 0, add = "jitter")
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_bray_season_within_bysite.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# ---- Within-site distances, broken down by season ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "site", measure = "bray")
t1$cal_group_distance(within_group = TRUE, by_group = "season")
write.csv(t1$res_group_distance,
	file.path(output_dir, "BetaDiv_rarefy_bray_site_within_byseason.csv"))
t1$cal_group_distance_diff(method = "wilcox")
write.csv(t1$res_group_distance_diff,
	file.path(output_dir, "BetaDiv_rarefy_bray_site_within_byseason_diff.csv"))
g1 <- t1$plot_group_distance(xtext_angle = 0, add = "jitter")
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_bray_site_within_byseason.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# ---- Between-season distances, broken down by site ----
# Note: with only 2 seasons there is 1 between-group pair (February_September),
# so cal_group_distance_diff has nothing to compare and is skipped here.
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "bray")
t1$cal_group_distance(within_group = FALSE, by_group = "site")
g1 <- t1$plot_group_distance(xtext_angle = 0, add = "jitter") +
	theme(legend.title = element_blank(), legend.text = element_text(size = 12))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_bray_season_between_bysite.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)

# ---- Between-site distances, broken down by season ----
# Note: with only 2 sites there is 1 between-group pair (ColaDeBallena_EsteroZacatecas),
# so cal_group_distance_diff has nothing to compare and is skipped here.
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "site", measure = "bray")
t1$cal_group_distance(within_group = FALSE, by_group = "season")
g1 <- t1$plot_group_distance(xtext_angle = 0, add = "jitter") +
	theme(legend.title = element_blank(), legend.text = element_text(size = 12))
cowplot::save_plot(file.path(output_dir, "BetaDiv_rarefy_bray_site_between_byseason.png"),
	g1, base_aspect_ratio = 1.25, dpi = 300, base_height = 6)


###########################################################################################
# perMANOVA
###########################################################################################

# ---- one-way: season ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_manova(manova_all = TRUE, group = "season")
write.csv(t1$res_manova,
	file.path(output_dir, "BetaDiv_rarefy_bray_perMANOVA_season.csv"))

# ---- one-way: site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_manova(manova_all = TRUE, group = "site")
write.csv(t1$res_manova,
	file.path(output_dir, "BetaDiv_rarefy_bray_perMANOVA_site.csv"))

# ---- one-way: combined group (pairwise) ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_manova(manova_all = FALSE, group = "group")
write.csv(t1$res_manova,
	file.path(output_dir, "BetaDiv_rarefy_bray_perMANOVA_group_pairwise.csv"))

# ---- two-way: season * site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_manova(manova_set = "season*site")
write.csv(t1$res_manova,
	file.path(output_dir, "BetaDiv_rarefy_bray_perMANOVA_twoway.csv"))

# ---- one-way: season, tested separately within each site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_manova(manova_all = FALSE, group = "season", by_group = "site")
write.csv(t1$res_manova,
	file.path(output_dir, "BetaDiv_rarefy_bray_perMANOVA_season_bysite.csv"))


###########################################################################################
# ANOSIM
###########################################################################################

# ---- overall: season ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_anosim(paired = FALSE, group = "season")
write.csv(t1$res_anosim,
	file.path(output_dir, "BetaDiv_rarefy_bray_ANOSIM_season.csv"))

# ---- overall: site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_anosim(paired = FALSE, group = "site")
write.csv(t1$res_anosim,
	file.path(output_dir, "BetaDiv_rarefy_bray_ANOSIM_site.csv"))

# ---- pairwise: combined group ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_anosim(paired = TRUE, group = "group")
write.csv(t1$res_anosim,
	file.path(output_dir, "BetaDiv_rarefy_bray_ANOSIM_group_pairwise.csv"))

# ---- season within each site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, measure = "bray")
t1$cal_anosim(paired = TRUE, group = "season", by_group = "site")
write.csv(t1$res_anosim,
	file.path(output_dir, "BetaDiv_rarefy_bray_ANOSIM_season_bysite.csv"))


###########################################################################################
# PERMDISP — test homogeneity of group dispersions
###########################################################################################

# ---- by season ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "season", measure = "bray")
t1$cal_betadisper()
capture.output(t1$res_betadisper,
	file = file.path(output_dir, "BetaDiv_rarefy_bray_permdisp_season.txt"))

# ---- by site ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "site", measure = "bray")
t1$cal_betadisper()
capture.output(t1$res_betadisper,
	file = file.path(output_dir, "BetaDiv_rarefy_bray_permdisp_site.txt"))

# ---- by combined group ----
t1 <- trans_beta$new(dataset = tmp_rarefy, group = "group", measure = "bray")
t1$cal_betadisper()
capture.output(t1$res_betadisper,
	file = file.path(output_dir, "BetaDiv_rarefy_bray_permdisp_group.txt"))
