##
## Correlation analysis between differentially abundant families and environmental factors
## Groups: season (February / September) and site (ColaDeBallena / EsteroZacatecas)
##
## Water variables:    water_COND, water_TDS, water_ORP, water_pH, water_salinity, water_temperature
## # evaluated but did decide NOT to use them*: Sediment variables: sed_pH, sed_EC, sed_OM, sed_orgN, sed_P, sed_K, sed_sand, sed_silt, sed_clay
## #*sediment results were mainly due to site differences since no seasonal data 
## # was available. Consequently, excluded from beta diversity analysis.


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
library(readxl)
theme_set(theme_bw())

######################################################
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}

# load preprocessed microtable (family-filtered, consistent with Steps 12-14)
load(file.path(output_dir, "tmp_microtable_preproc.RData"))

# refresh sample_table with full environmental data from Sample_info.xlsx
# (the preproc object was saved before water/sediment columns were added)
tmp_sample <- as.data.frame(
	read_excel("Input/1.Amplicon/QIIME2_qza/Sample_info.xlsx"),
	stringsAsFactors = FALSE
)
rownames(tmp_sample) <- tmp_sample[[1]]
tmp_sample <- tmp_sample[, -1]
tmp_microtable$sample_table <- tmp_sample[rownames(tmp_microtable$sample_table), ]
tmp_microtable$tidy_dataset()

# Column indices in the refreshed sample_table:
# season(1) site(2) tree.id(3) transect.position(4) species(5) sample.type(6)
# latitude(7) longitude(8) collection.date(9)
# water_COND(10) water_TDS(11) water_ORP(12) water_pH(13) water_salinity(14) water_temperature(15)
# sed_pH(16) sed_EC(17) sed_OM(18) sed_orgN(19) sed_P(20) sed_K(21) sed_sand(22) sed_silt(23) sed_clay(24)

water_cols   <- 10:15
sed_cols     <- 16:24
all_env_cols <- 10:24

# calculate relative abundance
tmp_microtable$cal_abund(rel = TRUE)


######################################################
# Step 1: identify differentially abundant families
# using glmm_beta with season + site formula
######################################################

formula  <- "season+site"
taxlevel <- "Family"

tmp_transdiff <- trans_diff$new(dataset = tmp_microtable, method = "glmm_beta",
	formula = formula, taxa_level = taxlevel, filter_thres = 0.001)

# select families with at least one significant factor (excluding intercept)
select_taxa <- tmp_transdiff$res_diff %>%
	.[!is.na(.$Estimate), ] %>%
	.[.$Factors != "(Intercept)", ] %>%
	.[.$Factors != "Model", ] %>%
	.[grepl("*", .$Significance, fixed = TRUE), ] %>%
	.$Taxa %>%
	unique

message(length(select_taxa), " families selected for environmental correlation.")


######################################################
# Step 2: correlate significant families with
#         environmental variables
######################################################

## 2a. Water variables ---------------------------------------------------
tmp_transenv <- trans_env$new(dataset = tmp_microtable,
	env_cols = water_cols, standardize = TRUE)
tmp_transenv$cal_cor(use_data = "other", p_adjust_method = "fdr",
	other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor,
	file.path(output_dir, "Diff_abund_env_cor_water.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_water_cluster.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)


## 2b. Sediment variables ------------------------------------------------
tmp_transenv <- trans_env$new(dataset = tmp_microtable,
	env_cols = sed_cols, standardize = TRUE)
tmp_transenv$cal_cor(use_data = "other", p_adjust_method = "fdr",
	other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor,
	file.path(output_dir, "Diff_abund_env_cor_sediment.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_sediment_cluster.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)


## 2c. All environmental variables together ------------------------------
tmp_transenv <- trans_env$new(dataset = tmp_microtable,
	env_cols = all_env_cols, standardize = TRUE)
tmp_transenv$cal_cor(use_data = "other", p_adjust_method = "fdr",
	other_taxa = select_taxa)
write.csv(tmp_transenv$res_cor,
	file.path(output_dir, "Diff_abund_env_cor_allenv.csv"))

g1 <- tmp_transenv$plot_cor(cluster_ggplot = "both")
cowplot::save_plot(file.path(output_dir, "Diff_abund_env_cor_allenv_cluster.png"),
	g1, base_aspect_ratio = 1.2, dpi = 300, base_height = 8)


######################################################
# Step 3: visualize the differential test result
#         alongside the correlation heatmap
######################################################

# filter res_diff to selected taxa and meaningful factors only
tmp_res <- tmp_transdiff$res_diff %>%
	.[.$Taxa %in% select_taxa, ] %>%
	.[!.$Factors %in% c("(Intercept)", "Model"), ]

# clean factor and taxa labels for display
tmp_res$Factors %<>% gsub("seasonSeptember",     "Season: September",      .)
tmp_res$Factors %<>% gsub("siteEsteroZacatecas", "Site: Estero Zacatecas", .)
tmp_res$Taxa    %<>% gsub(".*f__", "", .)

tmp_transdiff$res_diff <- tmp_res
write.csv(tmp_transdiff$res_diff,
	file.path(output_dir, "Diff_abund_Family_glmmbeta_significant.csv"))

g1 <- tmp_transdiff$plot_diff_bar(
	heatmap_cell    = "Estimate",
	heatmap_lab_fill = "Betareg\nCoef",
	cluster_ggplot  = "both",
	xtext_angle     = 30)
cowplot::save_plot(file.path(output_dir, "Diff_abund_Family_glmmbeta_heatmap.png"),
	g1, base_aspect_ratio = 1.0, dpi = 300, base_height = 8)
