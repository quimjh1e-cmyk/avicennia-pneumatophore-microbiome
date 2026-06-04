##
## use trans_venn class to analyze shared and unique core ASVs among groups (season x site)
##
## Taxonomic composition of intersections is reported at both Family and Genus level.
## Family is the primary level (69.7% ASV coverage vs 36.6% for Genus).
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
######################################################
# check output directory
output_dir <- "./Output/1.Amplicon/Stage3_coremicrobiome"
if(! dir.exists(output_dir)){
	stop("Please first run the script in last step !")
}
load(file.path(output_dir, "Coretaxa_February_ColaDeBallena.RData"))
load(file.path(output_dir, "Coretaxa_February_EsteroZacatecas.RData"))
load(file.path(output_dir, "Coretaxa_September_ColaDeBallena.RData"))
load(file.path(output_dir, "Coretaxa_September_EsteroZacatecas.RData"))

# load rarefied microtable object
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
load(input_path)
######################################################

# generate a total core-taxa microtable object (union of all groups)
tmp_microtable <- clone(amplicon_16S_microtable_rarefy)
tmp_microtable$otu_table %<>% .[rownames(.) %in% c(FCBA$taxa_names(), FEZA$taxa_names(), SCBA$taxa_names(), SEZA$taxa_names()), ]
tmp_microtable$tidy_dataset()
# merge samples into one per group using the season+site combination
# create a combined grouping column
tmp_microtable$sample_table$group <- paste(tmp_microtable$sample_table$season, tmp_microtable$sample_table$site, sep = "_")
tmp_merge <- tmp_microtable$merge_samples("group")

# analyze intersection set with trans_venn class
trans_vennobj <- trans_venn$new(tmp_merge, ratio = "numratio", name_joint = "-")
# remove the groups with 0
trans_vennobj$data_summary %<>% .[.$Counts != 0, ]
# output data_summary to file
# The "Count" column represents the feature number in each set; The "Abundance" column denotes the percentage that each part constitutes of the total.
write.csv(trans_vennobj$data_summary, file.path(output_dir, "Coretaxa_intersection.csv"))

# Season colour scheme  — blue = February, red = September
# Site colour scheme    — light orange = Cola de Ballena, light green = Estero Zacatecas
# Combined (group_colors): season sets the hue family, site sets light vs. dark within it
group_colors <- c(
	"February_ColaDeBallena"    = "#444444",   # light blue  (Feb + CB)
	"February_EsteroZacatecas"  = "#444444",   # dark blue   (Feb + EZ)
	"September_ColaDeBallena"   = "#444444",   # light red   (Sep + CB)
	"September_EsteroZacatecas" = "#444444"    # dark red    (Sep + EZ)
)
# Desaturated versions used as row-background stripes (keep dots readable)
group_bg_colors <- c(
	"February_ColaDeBallena"    = "#f4a070",   # pale blue
	"February_EsteroZacatecas"  = "#70b870",   # medium blue
	"September_ColaDeBallena"   = "#f4a070",   # pale red
	"September_EsteroZacatecas" = "#70b870"    # medium red
)
# Site-only colour scheme (for site-level plots throughout the project)
site_colors <- c(
	"ColaDeBallena"   = "#f4a070",   # light orange
	"EsteroZacatecas" = "#70b870"    # light green
)

# Figure 2a — UpSet-style intersection bar chart
# left_bar_fill   : colours the per-group set-size bars (left panel)
# bottom/left _background_fill : soft row stripes so each group is immediately
#                                identifiable in the dot matrix and left panel
# up_bar_fill     : single neutral colour for intersection-size bars (no per-
#                   intersection group mapping is available in this API)
g1 <- trans_vennobj$plot_bar(
	sort_samples             = FALSE,
	up_y_title               = "Intersection size",
	up_y_title_size          = 12,
	up_y_text_size           = 9,
	up_bar_fill              = "#444444",
	up_bar_width             = 0.65,
	bottom_height            = 0.9,
	bottom_point_size        = 5,
	bottom_point_color       = "#1a1a1a",
	bottom_y_text_size       = 11,
	bottom_line_width        = 1.0,
	bottom_line_colour       = "#1a1a1a",
	bottom_background_fill   = group_bg_colors,
	bottom_background_alpha  = 0.9,
	left_width               = 0.38,
	left_bar_fill            = group_colors,
	left_bar_width           = 0.65,
	left_bar_alpha           = 0.9,
	left_x_text_size         = 9,
	left_background_fill     = group_bg_colors,
	left_background_alpha    = 0.9
)
cowplot::save_plot(file.path(output_dir, "Coretaxa_intersection.png"),
	g1, base_aspect_ratio = 1.4, dpi = 300, base_height = 7)


# transform intersection data to microtable object for other visualization
trans_vennobj_mt <- trans_vennobj$trans_comm(use_frequency = TRUE)
# calculate relative abundance at all taxonomic levels
trans_vennobj_mt$cal_abund()

# Ordered elements in x-axis according to count — keeps the order consistent with
# 'Coretaxa_intersection.png' in the composition bar plots below
tmp_x <- trans_vennobj$data_summary %>%
	tibble::rownames_to_column(.) %>%
	.[order(.$Counts, decreasing = TRUE), ] %>%
	.$rowname

# Helper: strip full lineage path from taxon labels in a trans_abund object.
# microeco may store taxonomy as "p__X|c__Y|f__FamilyName" or with rank prefixes.
# The taxonomy column may be named by taxrank or differently — detect it robustly.
strip_taxon_labels <- function(abund_obj) {
	df  <- abund_obj$data_abund
	# Identify the taxonomy column: character columns that are not numeric
	# and not the sample-ID column (which typically contains our sample names)
	sample_ids <- rownames(abund_obj$dataset$sample_table)
	char_cols  <- names(df)[sapply(df, is.character)]
	tax_col    <- char_cols[!sapply(char_cols, function(cn) any(df[[cn]] %in% sample_ids))]
	if (length(tax_col) == 0) tax_col <- char_cols   # fallback: use all char cols
	tax_col <- tax_col[1]
	df[[tax_col]] <- gsub("^.*\\|", "", df[[tax_col]])   # strip lineage prefix
	df[[tax_col]] <- gsub("^[a-z]__", "", df[[tax_col]]) # strip rank prefix (f__, g__, …)
	abund_obj$data_abund <- df
	abund_obj
}

# Loop over Family (primary: 69.7 % ASV coverage) and Genus (secondary: 36.6 %)
for (taxrank in c("Family", "Genus")) {

	# CSV: proportion of core ASVs within each intersection set attributed to each taxon
	# Each column is an intersection group; each row is a taxon at the chosen rank.
	write.csv(
		trans_vennobj_mt$taxa_abund[[taxrank]],
		file.path(output_dir, paste0("Coretaxa_intersection_", taxrank, "ratio.csv"))
	)

	# Bar plot: top-10 taxa at this rank across intersection sets
	trans_vennobj_mt_abund <- trans_abund$new(
		dataset = trans_vennobj_mt, taxrank = taxrank, ntaxa = 10)
	trans_vennobj_mt_abund <- strip_taxon_labels(trans_vennobj_mt_abund)

	g_tax <- trans_vennobj_mt_abund$plot_bar(
		bar_full          = TRUE,
		legend_text_italic = (taxrank == "Genus"),   # italics only for genus
		xtext_angle       = 40,
		order_x           = tmp_x)
	g_tax <- g_tax +
		ylab("Ratio (%)") +
		ggtitle(paste0("Core taxa intersections — ", taxrank, " level")) +
		theme(plot.margin = unit(c(1, 0, 0, 2), "cm"))

	cowplot::save_plot(
		file.path(output_dir, paste0("Coretaxa_intersection_", taxrank, "ratio.png")),
		g_tax, base_aspect_ratio = 1.3, dpi = 300, base_height = 7)
}


##########################################################################
# Per-sample stacked bar: top 10 families across all 24 samples
# Uses the full rarefied community (not filtered to core taxa).
##########################################################################

# Order samples: February before September, Cola de Ballena before Estero Zacatecas
tmp_full <- clone(amplicon_16S_microtable_rarefy)
tmp_full$sample_table <- tmp_full$sample_table[
	order(tmp_full$sample_table$season, tmp_full$sample_table$site), ]
tmp_full$tidy_dataset()
tmp_full$cal_abund()

abund_fam <- trans_abund$new(dataset = tmp_full, taxrank = "Genus", ntaxa = 20)
abund_fam <- strip_taxon_labels(abund_fam)

g_samples <- abund_fam$plot_bar(
	bar_full   = TRUE,
	xtext_angle = 45,
	legend_text_italic = FALSE)
g_samples <- g_samples +
	ylab("Relative abundance (%)") +
	xlab("") +
	ggtitle("Top 20 genera — all samples") +
	theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))

cowplot::save_plot(
	file.path(output_dir, "Top20Genera_perSample.png"),
	g_samples, base_aspect_ratio = 2.0, dpi = 300, base_height = 6)
