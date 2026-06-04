##
## Correlation between alpha diversity and environmental variables
##
## Water variables:   water_COND, water_TDS, water_ORP, water_pH, water_salinity, water_temperature
## # evaluated but did decide NOT to use them*: Sediment variables: sed_pH, sed_EC, sed_OM, sed_orgN, sed_P, sed_K, sed_sand, sed_silt, sed_clay
##
##   (sediment measured in September only; assumed stable and applied to both seasons)
##
## #*sediment results were mainly due to site differences since no seasonal data 
## # was available. Consequently, excluded from beta diversity analysis.
##
## Groups: season (February, September) x site (ColaDeBallena, EsteroZacatecas)
## With n = 6 per group, Spearman correlation is used (non-parametric, robust for small n).
##


###########################
# load packages
library(microeco)
library(magrittr)
library(readxl)
###########################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage4_AlphaDiversity"

# load data
input_path <- file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData")
if(! file.exists(input_path)){
	stop("Please first run the script in step7 !")
}
load(input_path)
###########################

# Refresh sample_table from updated Sample_info.xlsx (includes water_ and sed_ columns).
# This is necessary because the RData file was saved before the sediment data was added.
tmp_sample <- as.data.frame(
	read_excel("Input/1.Amplicon/QIIME2_qza/Sample_info.xlsx"),
	stringsAsFactors = FALSE
)
rownames(tmp_sample) <- tmp_sample[[1]]
tmp_sample <- tmp_sample[, -1]   # drop sample-id column (now stored as rownames only)
tmp_microtable$sample_table <- tmp_sample[rownames(tmp_microtable$sample_table), ]

# Column indices in the refreshed sample_table:
# season(1) site(2) tree-id(3) transect-position(4) species(5) sample-type(6)
# latitude(7) longitude(8) collection-date(9)
# water_COND(10) water_TDS(11) water_ORP(12) water_pH(13) water_salinity(14) water_temperature(15)
# sed_pH(16) sed_EC(17) sed_OM(18) sed_orgN(19) sed_P(20) sed_K(21) sed_sand(22) sed_silt(23) sed_clay(24)

water_cols   <- 10:15   # intertidal water variables
sed_cols     <- 16:24   # sediment variables
all_env_cols <- 10:24   # all environmental variables combined
spatial_cols <- 7:8     # latitude, longitude



# Helper: build trans_env, run Spearman correlation, write CSV and save heatmap
run_alphacor <- function(dataset, env_cols, label, height = 6, aspect = 1.5) {
	env_obj <- trans_env$new(dataset = dataset, env_cols = env_cols)
	env_obj$cal_cor(add_abund_table = dataset$alpha_diversity, cor_method = "spearman")
	write.csv(env_obj$res_cor,
		file.path(output_dir, paste0("AlphaDiv_Env_", label, "_spearman.csv")))
	g <- env_obj$plot_cor()
	cowplot::save_plot(
		file.path(output_dir, paste0("AlphaDiv_Env_", label, "_spearman_heatmap.png")),
		g, base_aspect_ratio = aspect, dpi = 300, base_height = height)
}
# Strip water_ prefix for display — restore immediately after
#colnames(tmp_microtable$sample_table)[water_cols] <- sub("^water_", "",
                                                         #colnames(tmp_microtable$sample_table)[water_cols])
#run_alphacor(tmp_microtable, water_cols, "allsamples_water")
#colnames(tmp_microtable$sample_table)[water_cols] <- paste0("water_",
                                                            #colnames(tmp_microtable$sample_table)[water_cols])

##########################################################################
# 1. All samples (n = 24)
##########################################################################

run_alphacor(tmp_microtable, water_cols,   "allsamples_water")
run_alphacor(tmp_microtable, sed_cols,     "allsamples_sediment",  aspect = 1.5)
run_alphacor(tmp_microtable, all_env_cols, "allsamples_allenv",    aspect = 1.8)


##########################################################################
# 2. By season — February (n = 12) and September (n = 12)
##########################################################################

for(seas in c("February", "September")){
	tmp_sub <- clone(tmp_microtable)
	tmp_sub$sample_table %<>% .[.$season == seas, ]
	tmp_sub$tidy_dataset()

	run_alphacor(tmp_sub, water_cols, paste0(seas, "_water"))
	run_alphacor(tmp_sub, sed_cols,   paste0(seas, "_sediment"), aspect = 1.5)
}


##########################################################################
# 3. By site — ColaDeBallena (n = 12) and EsteroZacatecas (n = 12)
##########################################################################

for(st in c("ColaDeBallena", "EsteroZacatecas")){
	tmp_sub <- clone(tmp_microtable)
	tmp_sub$sample_table %<>% .[.$site == st, ]
	tmp_sub$tidy_dataset()

	run_alphacor(tmp_sub, water_cols, paste0(st, "_water"))
	run_alphacor(tmp_sub, sed_cols,   paste0(st, "_sediment"), aspect = 1.5)
}


##########################################################################
# 4. By group (season × site, n = 6 each)
##########################################################################

groups <- list(
	February_ColaDeBallena    = list(season = "February",  site = "ColaDeBallena"),
	February_EsteroZacatecas  = list(season = "February",  site = "EsteroZacatecas"),
	September_ColaDeBallena   = list(season = "September", site = "ColaDeBallena"),
	September_EsteroZacatecas = list(season = "September", site = "EsteroZacatecas")
)

for(grp_name in names(groups)){
	seas <- groups[[grp_name]]$season
	st   <- groups[[grp_name]]$site

	tmp_sub <- clone(tmp_microtable)
	tmp_sub$sample_table %<>% .[.$season == seas & .$site == st, ]
	tmp_sub$tidy_dataset()

	run_alphacor(tmp_sub, water_cols, paste0(grp_name, "_water"))
	run_alphacor(tmp_sub, sed_cols,   paste0(grp_name, "_sediment"), aspect = 1.5)
}


##########################################################################
# 5. Spatial correlation — all samples (n = 24), latitude + longitude
##########################################################################

tmp <- trans_env$new(dataset = tmp_microtable, env_cols = spatial_cols)
tmp$cal_cor(
	add_abund_table = tmp_microtable$alpha_diversity,
	cor_method = "spearman"
)
write.csv(tmp$res_cor, file.path(output_dir, "AlphaDiv_Spatial_allsamples_spearman.csv"))

g1 <- tmp$plot_cor()
cowplot::save_plot(
	file.path(output_dir, "AlphaDiv_Spatial_allsamples_spearman_heatmap.png"),
	g1, base_aspect_ratio = 1.0, dpi = 300, base_height = 6)
