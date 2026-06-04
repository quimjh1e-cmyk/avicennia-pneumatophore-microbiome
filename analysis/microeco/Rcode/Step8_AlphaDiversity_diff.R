##
## Differential test for alpha diversity
##


###########################
# load packages
library(microeco)
library(magrittr)
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage4_AlphaDiversity"

# load data
input_path <- file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData")
if(! file.exists(input_path)){
	stop("Please first run the script of last step !")
}
load(input_path)
###########################


measure <- "Shannon"


##########################################################################
# 1. Compare seasons (February vs September) — Wilcoxon test
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "season")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
# errorbar plot
g1 <- tmp$plot_alpha(measure = measure, plot_type = "errorbar", add_line = TRUE, line_type = 2)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_wilcox_errorbar.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

g1

##########################################################################
# 2. Compare sites (ColaDeBallena vs EsteroZacatecas) — Wilcoxon test
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "site")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_site_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_site_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


##########################################################################
# 3. Compare seasons within each site — Wilcoxon test (by_group)
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "season", by_group = "site")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_bysite_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_bysite_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


##########################################################################
# 4. Compare sites within each season — Wilcoxon test (by_group)
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "site", by_group = "season")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_site_byseason_wilcox.csv"))
g1 <- tmp$plot_alpha(measure = measure)
cowplot::save_plot(file.path(output_dir, "AlphaDiv_site_byseason_wilcox_boxplot.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)


##########################################################################
# 5. Two-way ANOVA: season * site (interaction + main effects)
##########################################################################

# for Shannon only
tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "anova", measure = measure, formula = "season*site")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_site_twowayanova.csv"))

# for all alpha diversity measures
tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "anova", formula = "season*site")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_site_twowayanova_allmeasures.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_site_twowayanova_allmeasures_heatmap.png"), g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 7)


##########################################################################
# 6. Linear model: season + site
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable)
tmp$cal_diff(method = "lm", formula = "season+site")
# filter out the intercept
tmp$res_diff %<>% .[.$Factors != "(Intercept)", ]
# clean up factor labels for plotting
tmp$res_diff$Factors %<>% gsub("seasonSeptember", "Season: September", .)
tmp$res_diff$Factors %<>% gsub("siteEsteroZacatecas", "Site: EsteroZacatecas", .)

g1 <- tmp$plot_alpha(heatmap_cell = "Estimate", heatmap_lab_fill = "Coef", xtext_angle = 30,
	text_x_order = c("Season: September", "Site: EsteroZacatecas")) +
	theme(legend.position = "left")
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_site_lm_heatmap.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)

# errorbar for Shannon only
tmp$res_diff %<>% .[.$Measure == "Shannon", ]
g1 <- tmp$plot_alpha(coefplot_sig_pos = 3)
g1 <- g1 + theme(axis.text.y = element_text(size = 13, color = "black"), axis.title.x = element_text(size = 13, color = "black"))
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_site_lm_Shannon_errorbar.png"), g1, base_aspect_ratio = 1.6, dpi = 300, base_height = 7)


##########################################################################
# 7. Linear mixed-effects model: season + site + (1|transect-position)
#    transect-position treated as random effect to account for repeated
#    sampling at the same physical transect points across seasons/sites
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable)
# Note: column name uses dot notation as loaded into R
tmp$cal_diff(method = "lme", formula = "season+site+(1|`transect-position`)")
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_site_lme.csv"))
g1 <- tmp$plot_alpha()
cowplot::save_plot(file.path(output_dir, "AlphaDiv_season_site_lme.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)


##########################################################################
# 8. Paired Wilcoxon test: seasons paired by transect-position within site
#    Each transect position (1-6) is sampled in both seasons at each site
##########################################################################

tmp <- trans_alpha$new(dataset = tmp_microtable, group = "season", by_group = "site", by_ID = "transect-position")
tmp$cal_diff(method = "wilcox", measure = measure)
write.csv(tmp$res_diff, file.path(output_dir, "AlphaDiv_season_bysite_wilcox_paired.csv"))

