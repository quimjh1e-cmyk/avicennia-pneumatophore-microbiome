##
## Differential abundance test at ASV level
## Groups: season (February / September) and site (ColaDeBallena / EsteroZacatecas)
##


######################################################
# load packages
library(microeco)
library(magrittr)

######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# preprocess data
# All samples are composite pneumatophore samples — no compartment subsetting needed
tmp_microtable <- clone(amplicon_16S_microtable)
tmp_microtable$cal_abund(rel = TRUE)

# Remove ASVs with very low relative abundance
tmp_microtable$filter_taxa(rel_abund = 0.0001)

# Remove ASVs without meaningful family-level annotation
tmp_microtable$tax_table %<>% .[
  !is.na(.$Family) &
    !grepl("^f?__$", .$Family) &
    .$Family != "f__Incertae_Sedis",
]
tmp_microtable$tidy_dataset()

save(tmp_microtable, file = file.path(output_dir, "tmp_microtable_preproc.RData"), compress = TRUE)


#############################################
##  ASV level   single factor  — season    ##
##  (2 groups: February, September)        ##
#############################################

group    <- "season"
taxlevel <- "ASV"

# metastat (no object saved — not used in Step 13 intersection analysis)
method <- "metastat"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))

# metagenomeSeq
method <- "metagenomeSeq"
tmp_metagenomeSeq <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_metagenomeSeq$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_metagenomeSeq, file = file.path(output_dir, "ASV_season_metagenomeSeq.RData"), compress = TRUE)
g_vol <- tmp_metagenomeSeq$plot_volcano()
cowplot::save_plot(file.path(output_dir, "Diff_abund_ASV_season_metagenomeSeq_volcano.png"),
	g_vol, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g_ab <- tmp_metagenomeSeq$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_ASV_season_metagenomeSeq_abund.png"),
	g_ab, base_aspect_ratio = 1.3, dpi = 300, base_height = 8)

# ALDEx2_t
method <- "ALDEx2_t"
tmp_ALDEx2_t <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ALDEx2_t$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_ALDEx2_t, file = file.path(output_dir, "ASV_season_ALDEx2_t.RData"), compress = TRUE)

# ALDEx2_kw
method <- "ALDEx2_kw"
tmp_ALDEx2_kw <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ALDEx2_kw$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_ALDEx2_kw, file = file.path(output_dir, "ASV_season_ALDEx2_kw.RData"), compress = TRUE)

# DESeq2
method <- "DESeq2"
tmp_DESeq2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_DESeq2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_DESeq2, file = file.path(output_dir, "ASV_season_DESeq2.RData"), compress = TRUE)

# edgeR
method <- "edgeR"
tmp_edgeR <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_edgeR$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_edgeR, file = file.path(output_dir, "ASV_season_edgeR.RData"), compress = TRUE)

# ancombc2
method <- "ancombc2"
tmp_ancombc2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ancombc2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_ancombc2, file = file.path(output_dir, "ASV_season_ancombc2.RData"), compress = TRUE)

# linda
method <- "linda"
tmp_linda <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_linda$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_linda, file = file.path(output_dir, "ASV_season_linda.RData"), compress = TRUE)

# maaslin2
tmp_dir <- "tmp_maaslin2"
dir.create(tmp_dir, showWarnings = FALSE)
method <- "maaslin2"
tmp_maaslin2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel,
	standardize = FALSE, fixed_effects = group,
	tmp_input_maaslin2  = file.path(tmp_dir, "maaslin2_season_input"),
	tmp_output_maaslin2 = file.path(tmp_dir, "maaslin2_season_output"))
write.csv(tmp_maaslin2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, ".csv")))
save(tmp_maaslin2, file = file.path(output_dir, "ASV_season_maaslin2.RData"), compress = TRUE)

# GMPR + wilcox
tmp_norm <- trans_norm$new(tmp_microtable)
norm_obj <- tmp_norm$norm(method = "GMPR")
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)
method <- "wilcox"
tmp_GMPR_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_GMPR_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, "_GMPR.csv")))
save(tmp_GMPR_wilcox, file = file.path(output_dir, "ASV_season_wilcox_GMPR.RData"), compress = TRUE)

# Wrench + wilcox
tmp_norm <- trans_norm$new(tmp_microtable)
norm_obj <- tmp_norm$norm(method = "Wrench", condition = group)
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)
method <- "wilcox"
tmp_Wrench_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_Wrench_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_season_", method, "_Wrench.csv")))
save(tmp_Wrench_wilcox, file = file.path(output_dir, "ASV_season_wilcox_Wrench.RData"), compress = TRUE)


#############################################
##  ASV level   single factor  — site      ##
##  (2 groups: ColaDeBallena,              ##
##             EsteroZacatecas)            ##
#############################################

group    <- "site"
taxlevel <- "ASV"

# metastat (no object saved — not used in Step 13 intersection analysis)
method <- "metastat"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))

# metagenomeSeq
method <- "metagenomeSeq"
tmp_metagenomeSeq <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_metagenomeSeq$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_metagenomeSeq, file = file.path(output_dir, "ASV_site_metagenomeSeq.RData"), compress = TRUE)
g_vol <- tmp_metagenomeSeq$plot_volcano()
cowplot::save_plot(file.path(output_dir, "Diff_abund_ASV_site_metagenomeSeq_volcano.png"),
	g_vol, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g_ab <- tmp_metagenomeSeq$plot_diff_abund(add_sig = TRUE)
cowplot::save_plot(file.path(output_dir, "Diff_abund_ASV_site_metagenomeSeq_abund.png"),
	g_ab, base_aspect_ratio = 1.3, dpi = 300, base_height = 8)

# ALDEx2_t
method <- "ALDEx2_t"
tmp_ALDEx2_t <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ALDEx2_t$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_ALDEx2_t, file = file.path(output_dir, "ASV_site_ALDEx2_t.RData"), compress = TRUE)

# ALDEx2_kw
method <- "ALDEx2_kw"
tmp_ALDEx2_kw <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ALDEx2_kw$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_ALDEx2_kw, file = file.path(output_dir, "ASV_site_ALDEx2_kw.RData"), compress = TRUE)

# DESeq2
method <- "DESeq2"
tmp_DESeq2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_DESeq2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_DESeq2, file = file.path(output_dir, "ASV_site_DESeq2.RData"), compress = TRUE)

# edgeR
method <- "edgeR"
tmp_edgeR <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_edgeR$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_edgeR, file = file.path(output_dir, "ASV_site_edgeR.RData"), compress = TRUE)

# ancombc2
method <- "ancombc2"
tmp_ancombc2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_ancombc2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_ancombc2, file = file.path(output_dir, "ASV_site_ancombc2.RData"), compress = TRUE)

# linda
method <- "linda"
tmp_linda <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_linda$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_linda, file = file.path(output_dir, "ASV_site_linda.RData"), compress = TRUE)

# maaslin2
tmp_dir <- "tmp_maaslin2"
dir.create(tmp_dir, showWarnings = FALSE)
method <- "maaslin2"
tmp_maaslin2 <- trans_diff$new(dataset = tmp_microtable, method = method, group = group, taxa_level = taxlevel,
	standardize = FALSE, fixed_effects = group,
	tmp_input_maaslin2  = file.path(tmp_dir, "maaslin2_site_input"),
	tmp_output_maaslin2 = file.path(tmp_dir, "maaslin2_site_output"))
write.csv(tmp_maaslin2$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, ".csv")))
save(tmp_maaslin2, file = file.path(output_dir, "ASV_site_maaslin2.RData"), compress = TRUE)

# GMPR + wilcox
tmp_norm <- trans_norm$new(tmp_microtable)
norm_obj <- tmp_norm$norm(method = "GMPR")
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)
method <- "wilcox"
tmp_GMPR_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_GMPR_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, "_GMPR.csv")))
save(tmp_GMPR_wilcox, file = file.path(output_dir, "ASV_site_wilcox_GMPR.RData"), compress = TRUE)

# Wrench + wilcox
tmp_norm <- trans_norm$new(tmp_microtable)
norm_obj <- tmp_norm$norm(method = "Wrench", condition = group)
norm_obj$add_rownames2taxonomy("ASV")
norm_obj$cal_abund(rel = FALSE)
method <- "wilcox"
tmp_Wrench_wilcox <- trans_diff$new(dataset = norm_obj, method = method, group = group, taxa_level = taxlevel)
write.csv(tmp_Wrench_wilcox$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_site_", method, "_Wrench.csv")))
save(tmp_Wrench_wilcox, file = file.path(output_dir, "ASV_site_wilcox_Wrench.RData"), compress = TRUE)


#############################################
##  ASV level   multiple factors           ##
##  season + site                          ##
#############################################

formula  <- "season+site"
taxlevel <- "ASV"

# DESeq2
method <- "DESeq2"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_multifactor_", method, ".csv")))

# ancombc2
method <- "ancombc2"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = NULL, fix_formula = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_multifactor_", method, ".csv")))

# linda
method <- "linda"
tmp <- trans_diff$new(dataset = tmp_microtable, method = method, group = formula, taxa_level = taxlevel)
write.csv(tmp$res_diff, file.path(output_dir, paste0("Diff_abund_ASV_multifactor_", method, ".csv")))
