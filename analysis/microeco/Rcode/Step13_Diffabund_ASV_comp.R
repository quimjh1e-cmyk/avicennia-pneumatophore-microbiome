##
## Differential abundance — method intersection and taxonomic composition at ASV level
## Two comparisons: season (February vs September) and site (ColaDeBallena vs EsteroZacatecas)
## Consensus-significant ASVs: present in >= 3 independent methods (adjusted p < 0.05)
## Taxonomic composition of intersection sets visualized at Family level
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
library(tibble)
theme_set(theme_bw())

######################################################
# output directory must already exist from Step 12
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	stop("Please first run Step 12!")
}

# load preprocessed microtable — reference for ASV names and taxonomy
load(file.path(output_dir, "tmp_microtable_preproc.RData"))
all_taxa_names <- rownames(tmp_microtable$otu_table)

######################################################
# Helper: build 0/1 presence-absence data.frame across methods
# diff_list : named list, one character vector of significant ASV names per method
# taxa_names: character vector of all ASV names (reference set)
build_pa_matrix <- function(diff_list, taxa_names) {
	tmp <- lapply(diff_list, function(x) {
		vec <- match(taxa_names, x)
		names(vec) <- taxa_names
		vec[is.na(vec)] <- 0
		vec[vec > 0]    <- 1
		vec
	})
	do.call(rbind, tmp) %>% t %>% as.data.frame
}


#############################################
##  Comparison 1: season                   ##
##  (February vs September)               ##
#############################################

season_files <- c(
	"ASV_season_metagenomeSeq.RData", "ASV_season_ALDEx2_t.RData",
	"ASV_season_ALDEx2_kw.RData",     "ASV_season_DESeq2.RData",
	"ASV_season_edgeR.RData",         "ASV_season_ancombc2.RData",
	"ASV_season_linda.RData",         "ASV_season_maaslin2.RData",
	"ASV_season_wilcox_GMPR.RData",   "ASV_season_wilcox_Wrench.RData"
)
missing <- season_files[!file.exists(file.path(output_dir, season_files))]
if(length(missing) > 0){
	stop("Missing season objects — re-run the season section of Step 12 first:\n",
	     paste(missing, collapse = "\n"))
}

load(file.path(output_dir, "ASV_season_metagenomeSeq.RData"))
load(file.path(output_dir, "ASV_season_ALDEx2_t.RData"))
load(file.path(output_dir, "ASV_season_ALDEx2_kw.RData"))
load(file.path(output_dir, "ASV_season_DESeq2.RData"))
load(file.path(output_dir, "ASV_season_edgeR.RData"))
load(file.path(output_dir, "ASV_season_ancombc2.RData"))
load(file.path(output_dir, "ASV_season_linda.RData"))
load(file.path(output_dir, "ASV_season_maaslin2.RData"))
load(file.path(output_dir, "ASV_season_wilcox_GMPR.RData"))
load(file.path(output_dir, "ASV_season_wilcox_Wrench.RData"))

diff_taxa_season <- list()
diff_taxa_season[["metagenomeSeq"]] <- tmp_metagenomeSeq$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["ALDEx2_t"]]      <- tmp_ALDEx2_t$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["ALDEx2_kw"]]     <- tmp_ALDEx2_kw$res_diff     %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["DESeq2"]]        <- tmp_DESeq2$res_diff        %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["edgeR"]]         <- tmp_edgeR$res_diff         %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["ancombc2"]]      <- tmp_ancombc2$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["linda"]]         <- tmp_linda$res_diff         %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["maaslin2"]]      <- tmp_maaslin2$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["GMPR_wilcox"]]   <- tmp_GMPR_wilcox$res_diff   %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_season[["Wrench_wilcox"]] <- tmp_Wrench_wilcox$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)

# build presence-absence matrix and microtable
pa_season <- build_pa_matrix(diff_taxa_season, all_taxa_names)
mt_season <- microtable$new(otu_table = pa_season, tax_table = tmp_microtable$tax_table)

# UpSet-style intersection bar plot
venn_season <- trans_venn$new(mt_season, ratio = "numratio", name_joint = "-")
venn_season$data_summary %<>% .[.$Counts > 2, ]   # show sets with > 2 ASVs

g1 <- venn_season$plot_bar(sort_samples = FALSE)
cowplot::save_plot(file.path(output_dir, "diff_ASV_season_methods_venn_bar.png"),
	g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 7)

# taxonomic composition of intersection sets at Family level
venn_season_mt <- venn_season$trans_comm(use_frequency = TRUE)
venn_season_mt$otu_table %<>% .[, colnames(.) %in% rownames(venn_season$data_summary)]
venn_season_mt$tidy_dataset()
venn_season_mt$cal_abund()
write.csv(venn_season_mt$taxa_abund$Family,
	file.path(output_dir, "diff_ASV_season_methods_venn_comp_Family.csv"))

tmp_x <- venn_season$data_summary %>%
	rownames_to_column %>%
	.[order(.$Counts, decreasing = TRUE), ] %>%
	.$rowname
venn_season_abund <- trans_abund$new(dataset = venn_season_mt, taxrank = "Family",
	ntaxa = 10, use_percentage = TRUE)
g2 <- venn_season_abund$plot_bar(bar_full = FALSE, legend_text_italic = TRUE,
		xtext_angle = 50, order_x = tmp_x) +
	ylab("Ratio (%)") +
	theme(legend.position = "left", plot.margin = unit(c(0, 0, 0, 4), "cm"))
cowplot::save_plot(file.path(output_dir, "diff_ASV_season_methods_venn_comp_Family.png"),
	g2, base_aspect_ratio = 1.5, dpi = 300, base_height = 8)

# remove season objects before loading site objects
rm(tmp_metagenomeSeq, tmp_ALDEx2_t, tmp_ALDEx2_kw, tmp_DESeq2, tmp_edgeR,
   tmp_ancombc2, tmp_linda, tmp_maaslin2, tmp_GMPR_wilcox, tmp_Wrench_wilcox)


#############################################
##  Comparison 2: site                     ##
##  (ColaDeBallena vs EsteroZacatecas)     ##
#############################################

load(file.path(output_dir, "ASV_site_metagenomeSeq.RData"))
load(file.path(output_dir, "ASV_site_ALDEx2_t.RData"))
load(file.path(output_dir, "ASV_site_ALDEx2_kw.RData"))
load(file.path(output_dir, "ASV_site_DESeq2.RData"))
load(file.path(output_dir, "ASV_site_edgeR.RData"))
load(file.path(output_dir, "ASV_site_ancombc2.RData"))
load(file.path(output_dir, "ASV_site_linda.RData"))
load(file.path(output_dir, "ASV_site_maaslin2.RData"))
load(file.path(output_dir, "ASV_site_wilcox_GMPR.RData"))
load(file.path(output_dir, "ASV_site_wilcox_Wrench.RData"))

diff_taxa_site <- list()
diff_taxa_site[["metagenomeSeq"]] <- tmp_metagenomeSeq$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["ALDEx2_t"]]      <- tmp_ALDEx2_t$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["ALDEx2_kw"]]     <- tmp_ALDEx2_kw$res_diff     %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["DESeq2"]]        <- tmp_DESeq2$res_diff        %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["edgeR"]]         <- tmp_edgeR$res_diff         %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["ancombc2"]]      <- tmp_ancombc2$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["linda"]]         <- tmp_linda$res_diff         %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["maaslin2"]]      <- tmp_maaslin2$res_diff      %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["GMPR_wilcox"]]   <- tmp_GMPR_wilcox$res_diff   %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)
diff_taxa_site[["Wrench_wilcox"]] <- tmp_Wrench_wilcox$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .$Taxa %>% gsub(".*\\|", "", .)

# build presence-absence matrix and microtable
pa_site <- build_pa_matrix(diff_taxa_site, all_taxa_names)
mt_site <- microtable$new(otu_table = pa_site, tax_table = tmp_microtable$tax_table)

# UpSet-style intersection bar plot
venn_site <- trans_venn$new(mt_site, ratio = "numratio", name_joint = "-")
venn_site$data_summary %<>% .[.$Counts > 2, ]

g1 <- venn_site$plot_bar(sort_samples = FALSE)
cowplot::save_plot(file.path(output_dir, "diff_ASV_site_methods_venn_bar.png"),
	g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 7)

# taxonomic composition of intersection sets at Family level
venn_site_mt <- venn_site$trans_comm(use_frequency = TRUE)
venn_site_mt$otu_table %<>% .[, colnames(.) %in% rownames(venn_site$data_summary)]
venn_site_mt$tidy_dataset()
venn_site_mt$cal_abund()
write.csv(venn_site_mt$taxa_abund$Family,
	file.path(output_dir, "diff_ASV_site_methods_venn_comp_Family.csv"))

tmp_x <- venn_site$data_summary %>%
	rownames_to_column %>%
	.[order(.$Counts, decreasing = TRUE), ] %>%
	.$rowname
venn_site_abund <- trans_abund$new(dataset = venn_site_mt, taxrank = "Family",
	ntaxa = 10, use_percentage = TRUE)
g2 <- venn_site_abund$plot_bar(bar_full = FALSE, legend_text_italic = TRUE,
		xtext_angle = 50, order_x = tmp_x) +
	ylab("Ratio (%)") +
	theme(legend.position = "left", plot.margin = unit(c(0, 0, 0, 4), "cm"))
cowplot::save_plot(file.path(output_dir, "diff_ASV_site_methods_venn_comp_Family.png"),
	g2, base_aspect_ratio = 1.5, dpi = 300, base_height = 8)
