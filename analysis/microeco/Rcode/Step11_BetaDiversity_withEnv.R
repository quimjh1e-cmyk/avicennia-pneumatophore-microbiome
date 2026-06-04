##
## Beta diversity analysis — environmental drivers
##
## Water variables:    water_COND, water_TDS, water_ORP, water_pH, water_salinity, water_temperature
## # evaluated but did decide NOT to use them*: Sediment variables: sed_pH, sed_EC, sed_OM, sed_orgN, sed_P, sed_K, sed_sand, sed_silt, sed_clay

##   (sediment measured in September only; assumed stable and applied to both seasons)

## #*sediment results were mainly due to site differences since no seasonal data 
## # was available. Consequently, excluded from beta diversity analysis.

## Ordination methods: dbRDA (distance-based), RDA (Genus level), CCA (Genus level)
## Additional analyses: Mantel test, scatter plot (Bray-Curtis vs. key variable)
## Grouping: season (February / September), site (ColaDeBallena / EsteroZacatecas)
##


###########################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
library(readxl)
theme_set(theme_bw())
###########################
# create output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage5_BetaDiversity"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load rarefied microtable
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
###########################
length(amplicon_16S_microtable_rarefy$rep_fasta)
# Refresh sample_table from updated Sample_info.xlsx (includes water_ and sed_ columns).
# This is necessary because the RData file was saved before the sediment data was added.
tmp_sample <- as.data.frame(
	read_excel("Input/1.Amplicon/QIIME2_qza/Sample_info.xlsx"),
	stringsAsFactors = FALSE
)
rownames(tmp_sample) <- tmp_sample[[1]]
tmp_sample <- tmp_sample[, -1]   # drop sample-id column (now stored as rownames only)
amplicon_16S_microtable_rarefy$sample_table <- tmp_sample[rownames(amplicon_16S_microtable_rarefy$sample_table), ]
amplicon_16S_microtable_rarefy$tidy_dataset()

# Column indices in the refreshed sample_table:
# season(1) site(2) tree-id(3) transect-position(4) species(5) sample-type(6)
# latitude(7) longitude(8) collection-date(9)
# water_COND(10) water_TDS(11) water_ORP(12) water_pH(13) water_salinity(14) water_temperature(15)
# sed_pH(16) sed_EC(17) sed_OM(18) sed_orgN(19) sed_P(20) sed_K(21) sed_sand(22) sed_silt(23) sed_clay(24)

water_cols   <- 10:15
sed_cols     <- 16:24
all_env_cols <- 10:24
spatial_cols <- 7:8

# Work with all 24 composite-pneumatophore samples
tmp_microtable <- clone(amplicon_16S_microtable_rarefy)

# Calculate Bray-Curtis beta diversity (required for dbRDA and Mantel test)
tmp_microtable$cal_betadiv(method = "bray")


#########################################################################################
## 1. dbRDA — distance-based RDA (Bray-Curtis)
#########################################################################################

method <- "dbRDA"

## 1a. Water variables ----------------------------------------------------------------
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = water_cols, standardize = TRUE)
t1$cal_ordination(method = method, use_measure = "bray")
t1$cal_ordination_anova()
t1$cal_ordination_envfit()
capture.output(t1$res_ordination_R2,
	file = file.path(output_dir, "BetaDiv_dbRDA_water_rawoutput.txt"))
capture.output(t1$res_ordination_envfit,
	file = file.path(output_dir, "BetaDiv_dbRDA_water_rawoutput.txt"), append = TRUE)
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_dbRDA_water_termssig.csv"))
write.csv(t1$res_ordination_axis,  file.path(output_dir, "BetaDiv_dbRDA_water_axissig.csv"))
t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
write.csv(t1$res_ordination_trans$df_sites,
	file.path(output_dir, "BetaDiv_dbRDA_water_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows,
	file.path(output_dir, "BetaDiv_dbRDA_water_trans_arrow.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_water.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# with forward feature selection
# env_data must exist in the calling frame: microeco's feature_sel uses eval(parent.frame())
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = water_cols, standardize = TRUE)
env_data <- t1$data_env
tryCatch({
	t1$cal_ordination(method = method, use_measure = "bray", feature_sel = TRUE)
	t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
	g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
	cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_water_featuresel.png"),
		g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
}, error = function(e) message("dbRDA water feature_sel: ", conditionMessage(e), " — plot skipped."))
rm(env_data)


## 1b. Sediment variables -------------------------------------------------------------
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = sed_cols, standardize = TRUE)
t1$cal_ordination(method = method, use_measure = "bray")
t1$cal_ordination_anova()
t1$cal_ordination_envfit()
capture.output(t1$res_ordination_R2,
	file = file.path(output_dir, "BetaDiv_dbRDA_sediment_rawoutput.txt"))
capture.output(t1$res_ordination_envfit,
	file = file.path(output_dir, "BetaDiv_dbRDA_sediment_rawoutput.txt"), append = TRUE)
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_dbRDA_sediment_termssig.csv"))
write.csv(t1$res_ordination_axis,  file.path(output_dir, "BetaDiv_dbRDA_sediment_axissig.csv"))
t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
write.csv(t1$res_ordination_trans$df_sites,
	file.path(output_dir, "BetaDiv_dbRDA_sediment_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows,
	file.path(output_dir, "BetaDiv_dbRDA_sediment_trans_arrow.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_sediment.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# with forward feature selection
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = sed_cols, standardize = TRUE)
env_data <- t1$data_env
tryCatch({
	t1$cal_ordination(method = method, use_measure = "bray", feature_sel = TRUE)
	t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
	g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
	cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_sediment_featuresel.png"),
		g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
}, error = function(e) message("dbRDA sediment feature_sel: ", conditionMessage(e), " — plot skipped."))
rm(env_data)


## 1c. All environmental variables ----------------------------------------------------
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = all_env_cols, standardize = TRUE)
t1$cal_ordination(method = method, use_measure = "bray")
t1$cal_ordination_anova()
t1$cal_ordination_envfit()
capture.output(t1$res_ordination_R2,
	file = file.path(output_dir, "BetaDiv_dbRDA_allenv_rawoutput.txt"))
capture.output(t1$res_ordination_envfit,
	file = file.path(output_dir, "BetaDiv_dbRDA_allenv_rawoutput.txt"), append = TRUE)
write.csv(t1$res_ordination_terms, file.path(output_dir, "BetaDiv_dbRDA_allenv_termssig.csv"))
write.csv(t1$res_ordination_axis,  file.path(output_dir, "BetaDiv_dbRDA_allenv_axissig.csv"))
t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
write.csv(t1$res_ordination_trans$df_sites,
	file.path(output_dir, "BetaDiv_dbRDA_allenv_trans_sample.csv"))
write.csv(t1$res_ordination_trans$df_arrows,
	file.path(output_dir, "BetaDiv_dbRDA_allenv_trans_arrow.csv"))
g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_allenv.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

# with forward feature selection
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = all_env_cols, standardize = TRUE)
env_data <- t1$data_env
tryCatch({
	t1$cal_ordination(method = method, use_measure = "bray", feature_sel = TRUE)
	t1$trans_ordination(adjust_arrow_length = TRUE, min_perc_env = 0.2, max_perc_env = 1)
	g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
	cowplot::save_plot(file.path(output_dir, "BetaDiv_dbRDA_allenv_featuresel.png"),
		g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
}, error = function(e) message("dbRDA allenv feature_sel: ", conditionMessage(e), " — plot skipped."))
rm(env_data)


#########################################################################################
## 2 & 3. RDA and CCA at Family and Genus level
##
## Family level (69.7 % ASV classification coverage) is the primary level.
## Genus level  (36.6 % coverage) is retained as supplementary.
## Both are run identically; output filenames encode the method and level.
#########################################################################################

# Helper: run one RDA/CCA block for a given method, taxa level, env set, and label
run_taxordination <- function(method, taxa_level, env_cols, env_label) {
	lbl <- paste0(method, "_", taxa_level, "_", env_label)
	t1 <- trans_env$new(dataset = tmp_microtable, env_cols = env_cols, standardize = TRUE)
	t1$cal_ordination(method = method, taxa_level = taxa_level)
	t1$cal_ordination_anova()
	t1$cal_ordination_envfit()
	capture.output(t1$res_ordination_R2,
		file = file.path(output_dir, paste0("BetaDiv_", lbl, "_rawoutput.txt")))
	capture.output(t1$res_ordination_envfit,
		file = file.path(output_dir, paste0("BetaDiv_", lbl, "_rawoutput.txt")), append = TRUE)
	write.csv(t1$res_ordination_terms,
		file.path(output_dir, paste0("BetaDiv_", lbl, "_termssig.csv")))
	write.csv(t1$res_ordination_axis,
		file.path(output_dir, paste0("BetaDiv_", lbl, "_axissig.csv")))
	t1$trans_ordination(adjust_arrow_length = TRUE,
		min_perc_env = 0.2, max_perc_env = 1, min_perc_tax = 0.2, max_perc_tax = 1)
	write.csv(t1$res_ordination_trans$df_sites,
		file.path(output_dir, paste0("BetaDiv_", lbl, "_trans_sample.csv")))
	write.csv(t1$res_ordination_trans$df_arrows,
		file.path(output_dir, paste0("BetaDiv_", lbl, "_trans_arrowenv.csv")))
	write.csv(t1$res_ordination_trans$df_arrows_spe,
		file.path(output_dir, paste0("BetaDiv_", lbl, "_trans_arrowspe.csv")))
	g1 <- t1$plot_ordination(plot_color = "season", plot_shape = "site")
	cowplot::save_plot(file.path(output_dir, paste0("BetaDiv_", lbl, ".png")),
		g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
}

# Family level first (primary — higher classification coverage: 69.7 %)
# then Genus level (supplementary — 36.6 % coverage)
for(taxa_lvl in c("Family", "Genus")){
	for(m in c("RDA", "CCA")){
		run_taxordination(m, taxa_lvl, water_cols, "water")
		run_taxordination(m, taxa_lvl, sed_cols,   "sediment")
	}
}


##########################################################################
## 4. Mantel test
##########################################################################

## 4a. Water variables ----------------------------------------------------------------
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = water_cols, standardize = TRUE)
t1$cal_mantel(use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_water_bray.csv"))

t1$cal_mantel(partial_mantel = TRUE, use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_water_bray_partial.csv"))

t1$cal_mantel(by_group = "season", use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_water_bray_bySeason.csv"))

t1$cal_mantel(by_group = "season", partial_mantel = TRUE, use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_water_bray_bySeason_partial.csv"))

t1$cal_mantel(by_group = "site", use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_water_bray_bySite.csv"))


## 4b. Sediment variables -------------------------------------------------------------
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = sed_cols, standardize = TRUE)
t1$cal_mantel(use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_sediment_bray.csv"))

t1$cal_mantel(partial_mantel = TRUE, use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_sediment_bray_partial.csv"))

t1$cal_mantel(by_group = "season", use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_sediment_bray_bySeason.csv"))

t1$cal_mantel(by_group = "site", use_measure = "bray")
write.csv(t1$res_mantel, file.path(output_dir, "BetaDiv_mantel_sediment_bray_bySite.csv"))


##########################################################################
## 5. Scatter plot — Bray-Curtis dissimilarity vs. key environmental variable
##########################################################################

## 5a. Bray-Curtis vs. water_salinity, coloured by season
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = water_cols, standardize = TRUE)
g1 <- t1$plot_scatterfit(
	x = tmp_microtable$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)],
	y = "water_salinity",
	type = "cor", group = "season", group_order = c("February", "September"),
	point_size = 3, point_alpha = 0.6, line_se = TRUE, line_size = 1.5, shape_values = c(16, 17),
	y_axis_title = "Euclidean distance of water salinity",
	x_axis_title = "Bray-Curtis distance", size = 5
)
cowplot::save_plot(file.path(output_dir, "BetaDiv_scatterfit_waterSalinity_Bray_bySeason.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

## 5b. Bray-Curtis vs. water_temperature, coloured by season
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = water_cols, standardize = TRUE)
g1 <- t1$plot_scatterfit(
  x = tmp_microtable$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)],
  y = "water_temperature",
  type = "cor", group = "season", group_order = c("February", "September"),
  point_size = 3, point_alpha = 0.6, line_se = TRUE, line_size = 1.5, shape_values = c(16, 17),
  y_axis_title = "Euclidean distance of water temperature",
  x_axis_title = "Bray-Curtis distance", size = 5
)
cowplot::save_plot(file.path(output_dir, "BetaDiv_scatterfit_waterTemperature_Bray_bySeason.png"),
                   g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)

## 5c. Bray-Curtis vs. sed_OM, coloured by site
t1 <- trans_env$new(dataset = tmp_microtable, env_cols = sed_cols, standardize = TRUE)
g1 <- t1$plot_scatterfit(
	x = tmp_microtable$beta_diversity$bray[rownames(t1$data_env), rownames(t1$data_env)],
	y = "sed_OM",
	type = "cor", group = "site", group_order = c("ColaDeBallena", "EsteroZacatecas"),
	point_size = 3, point_alpha = 0.6, line_se = TRUE, line_size = 1.5, shape_values = c(16, 17),
	y_axis_title = "Euclidean distance of sediment OM (%)",
	x_axis_title = "Bray-Curtis distance", size = 5
)
cowplot::save_plot(file.path(output_dir, "BetaDiv_scatterfit_sedOM_Bray_bySite.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 6)
