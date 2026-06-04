##
## Differential abundance test at Family level (aggregated reads)
## Groups: season (February / September) and site (ColaDeBallena / EsteroZacatecas)
## Methods: LEfSe, Wilcoxon + AST, ANOVA + AST, KW_dunn + AST, glmm_beta
## Also: all taxonomic levels overview and random-effects models
##
## Note: these relative-abundance-based methods are complementary to the
## count-based methods in Step 12. Family-level aggregation increases power
## by pooling reads across ASVs within each family.
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())

######################################################
# output directory
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}

# load preprocessed microtable (family-filtered, consistent with Steps 12-13)
load(file.path(output_dir, "tmp_microtable_preproc.RData"))
# calculate relative abundance required by all methods in this script
tmp_microtable$cal_abund(rel = TRUE)


##################################################################
##  Family level    relative abundance    single factor        ###
##  season (February vs September)                             ###
##################################################################

group    <- "season"
taxlevel <- "Family"

# LEfSe
method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_season_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_season_lefse_LDAbar.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_season_lefse_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# Wilcoxon + arcsine transformation
method <- "wilcox"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_season_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_season_wilcox_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# ANOVA + arcsine transformation
method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_season_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_season_ANOVA_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# glmm_beta
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = group,
	taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_season_", method, ".csv")))


##################################################################
##  Family level    relative abundance    single factor        ###
##  site (ColaDeBallena vs EsteroZacatecas)                    ###
##################################################################

group    <- "site"
taxlevel <- "Family"

# LEfSe
method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_site_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_site_lefse_LDAbar.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_site_lefse_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# Wilcoxon + arcsine transformation
method <- "wilcox"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_site_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_site_wilcox_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# ANOVA + arcsine transformation
method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_site_", method, ".csv")))
g1 <- tmp$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_site_ANOVA_abundsig.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)

# glmm_beta
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = group,
	taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_site_", method, ".csv")))


##################################################################
##  Family level    relative abundance    multiple factor      ###
##  season + site                                              ###
##################################################################

formula  <- "season+site"
taxlevel <- "Family"

# ANOVA + arcsine transformation
method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_multifactor_", method, ".csv")))

# linear regression + arcsine transformation
method <- "lm"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_multifactor_", method, ".csv")))

# glmm_beta + heatmap visualization of significant families
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_Family_multifactor_", method, ".csv")))

# filter results for visualization
tmp$res_diff %<>% .[.$Factors != "(Intercept)", ]
tmp$res_diff %<>% .[.$Factors != "Model", ]

# keep only families with at least one ** significant factor
tmp_sel <- c()
for(x in unique(tmp$res_diff$Taxa)){
	tmp_table <- tmp$res_diff[tmp$res_diff$Taxa == x, ]
	tmp_table <- tmp_table[tmp_table$Factors %in% c("seasonSeptember", "siteEsteroZacatecas"), ]
	if(any(grepl("**", tmp_table$Significance, fixed = TRUE))){
		tmp_sel <- c(tmp_sel, x)
	}
}
if(length(tmp_sel) > 0){
	tmp$res_diff %<>% .[tmp$res_diff$Taxa %in% tmp_sel, ]
	# clean labels
	tmp$res_diff$Factors %<>% gsub("seasonSeptember",     "Season: September",      .)
	tmp$res_diff$Factors %<>% gsub("siteEsteroZacatecas", "Site: Estero Zacatecas", .)
	tmp$res_diff$Taxa    %<>% gsub(".*f__", "", .)
	g2 <- tmp$plot_diff_bar(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef",
		xtext_angle = 30, xtext_size = 16, filter_feature = "",
		text_x_order = c("Season: September", "Site: Estero Zacatecas")) +
		theme(legend.position = "left")
	cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_multifactor_glmmbeta.png"),
		g2, base_aspect_ratio = 1.3, dpi = 300, base_height = 9)
} else {
	message("No families with ** significance in glmm_beta multifactor — heatmap skipped.")
}


##############################################################################
##  All taxonomic levels    relative abundance    single factor            ###
##  Provides an overview from phylum down to family                        ###
##############################################################################

taxlevel <- "all"

# --- season ---
group <- "season"

method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, alpha = 0.01, lefse_subgroup = NULL, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_season_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_allTax_season_lefse_LDAbar.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 5)
cowplot::save_plot(file.path(output_dir, "Diff_abund_allTax_season_lefse_cladogram.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 9)

method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_season_", method, ".csv")))

# --- site ---
group <- "site"

method <- "lefse"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, alpha = 0.01, lefse_subgroup = NULL, filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_site_", method, ".csv")))
g1 <- tmp$plot_diff_bar(width = 0.7)
cowplot::save_plot(file.path(output_dir, "Diff_abund_allTax_site_lefse_LDAbar.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
g1 <- tmp$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 5)
cowplot::save_plot(file.path(output_dir, "Diff_abund_allTax_site_lefse_cladogram.png"),
	g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 9)

method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_site_", method, ".csv")))


##############################################################################
##  All taxonomic levels    relative abundance    multiple factor          ###
##############################################################################

formula  <- "season+site"
taxlevel <- "all"

method <- "anova"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_multifactor_", method, ".csv")))

method <- "lm"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, transformation = "AST", filter_thres = 0.001)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_multifactor_", method, ".csv")))

method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, filter_thres = 0.005)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_allTax_multifactor_", method, ".csv")))


############## Random effects skipped due to software incompatibilities#############
############## Could not model tree-level variance #################################

##############################################################################
##  Random effects    Family level                                         ###
##  tree.id as random effect to account for tree-level variation           ###
##  Use if the same trees were sampled in both seasons (repeated measures) ###
##############################################################################

# rename column to avoid dot in formula parser
tmp_microtable$sample_table$tree_id <- tmp_microtable$sample_table[["tree.id"]]

formula  <- "season + site + (1|tree_id)"
taxlevel <- "Family"

# glmm_beta with random effect
# Conditional R2 = variance explained by fixed + random effects
# Marginal R2    = variance explained by fixed effects only
method <- "glmm_beta"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, formula = formula,
	taxa_level = taxlevel, filter_thres = 0.005)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_mixedeff_Family_", method, ".csv")))
