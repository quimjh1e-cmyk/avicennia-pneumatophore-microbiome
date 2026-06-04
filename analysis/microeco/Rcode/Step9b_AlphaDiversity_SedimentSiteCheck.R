##
## Diagnostic: are the February sediment–Simpson correlations driven by site separation?

## ## # evaluated but did decide NOT to use them*: Sediment variables: sed_pH, 
## ## # sed_EC, sed_OM, sed_orgN, sed_P, sed_K, sed_sand, sed_silt, sed_clay
## #*sediment results were mainly due to site differences since no seasonal data 
## # was available. Consequently, excluded from beta diversity analysis.

## The Spearman correlation (Step 9) tests whether, across the 12 February samples,
## a sediment gradient coincides with a Simpson gradient — but it does not know about
## site labels. This script checks whether CB and EZ have clearly different sed_P,
## sed_sand and sed_clay values (which would mean the correlation largely reflects a
## site effect) and visualises the Simpson gradient against those variables.
##
## Statistical test: Wilcoxon rank-sum (Mann-Whitney U), non-parametric, appropriate
## for n = 6 per group.
##


###########################
library(microeco)
library(magrittr)
library(readxl)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
theme_set(theme_bw())
###########################

output_dir <- "Output/1.Amplicon/Stage4_AlphaDiversity"

input_path <- file.path(output_dir, "amplicon_16S_microtable_rarefy_withalphadiv.RData")
if(!file.exists(input_path)) stop("Please first run the script in step7!")
load(input_path)

# Refresh sample_table
tmp_sample <- as.data.frame(
	read_excel("Input/1.Amplicon/QIIME2_qza/Sample_info.xlsx"),
	stringsAsFactors = FALSE
)
rownames(tmp_sample) <- tmp_sample[[1]]
tmp_sample <- tmp_sample[, -1]
tmp_microtable$sample_table <- tmp_sample[rownames(tmp_microtable$sample_table), ]

# Merge sample metadata with alpha diversity
meta      <- tmp_microtable$sample_table
alphadiv  <- tmp_microtable$alpha_diversity
df_all    <- cbind(meta, alphadiv[rownames(meta), ])

# Subset to February
df_feb <- df_all[df_all$season == "February", ]

# Sediment variables of interest (those significant in Step 9 February analysis)
sed_vars <- c("sed_P", "sed_sand", "sed_clay")

site_colors <- c("ColaDeBallena" = "#1f78b4", "EsteroZacatecas" = "#e66101")
site_order  <- c("ColaDeBallena", "EsteroZacatecas")
df_feb$site <- factor(df_feb$site, levels = site_order)


##########################################################################
## 1. Descriptive statistics and Wilcoxon tests — sed_P, sed_sand, sed_clay
##    between sites in February
##########################################################################

wilcox_results <- lapply(sed_vars, function(v) {
	cb <- df_feb[df_feb$site == "ColaDeBallena",  v]
	ez <- df_feb[df_feb$site == "EsteroZacatecas", v]
	wt <- wilcox.test(cb, ez, exact = FALSE)
	data.frame(
		variable   = v,
		CB_median  = median(cb, na.rm = TRUE),
		CB_IQR     = IQR(cb, na.rm = TRUE),
		EZ_median  = median(ez, na.rm = TRUE),
		EZ_IQR     = IQR(ez, na.rm = TRUE),
		W          = wt$statistic,
		p_value    = wt$p.value,
		stringsAsFactors = FALSE
	)
})
wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$p_adj <- p.adjust(wilcox_df$p_value, method = "BH")
write.csv(wilcox_df,
	file.path(output_dir, "SedimentSiteCheck_Wilcoxon_February.csv"),
	row.names = FALSE)


##########################################################################
## 2. Boxplots — sediment variables by site (February)
##########################################################################

# Pretty labels for plot axes
sed_labels <- c(
	sed_P    = "Phosphorus (ppm)",
	sed_sand = "Sand (%)",
	sed_clay = "Clay (%)"
)

# Build significance label from Wilcoxon p_adj
sig_label <- function(p) {
	if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else "ns"
}

box_plots <- lapply(sed_vars, function(v) {
	p_adj <- wilcox_df$p_adj[wilcox_df$variable == v]
	y_max <- max(df_feb[[v]], na.rm = TRUE)
	y_range <- diff(range(df_feb[[v]], na.rm = TRUE))

	ggplot(df_feb, aes(x = site, y = .data[[v]], fill = site)) +
		geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
		geom_jitter(width = 0.12, size = 2.5, alpha = 0.9, aes(color = site)) +
		annotate("text",
			x = 1.5, y = y_max + 0.1 * y_range,
			label = sig_label(p_adj), size = 6) +
		annotate("segment",
			x = 1, xend = 2,
			y = y_max + 0.05 * y_range,
			yend = y_max + 0.05 * y_range) +
		scale_fill_manual(values  = site_colors) +
		scale_color_manual(values = site_colors) +
		scale_x_discrete(labels = c("Cola de\nBallena", "Estero\nZacatecas")) +
		labs(y = sed_labels[v], x = NULL,
			title = sed_labels[v],
			caption = paste0("Wilcoxon p_adj = ", round(p_adj, 3))) +
		theme(legend.position = "none",
			plot.title = element_text(size = 11, face = "bold"),
			axis.text.x = element_text(size = 10))
})

g_boxes <- plot_grid(plotlist = box_plots, nrow = 1)
save_plot(
	file.path(output_dir, "SedimentSiteCheck_Boxplots_February.png"),
	g_boxes, ncol = 3, base_aspect_ratio = 0.9, dpi = 300, base_height = 5)


##########################################################################
## 3. Scatter plots — Simpson vs. sed_P, sed_sand, sed_clay (February)
##    Points coloured by site; smooth line across all 12 samples
##    This shows whether the Spearman correlation rides on the site split
##    or spans a genuine gradient within each site.
##########################################################################

scatter_plots <- lapply(sed_vars, function(v) {
	ggplot(df_feb, aes(x = .data[[v]], y = Simpson)) +
		# within-site trend lines
		geom_smooth(aes(color = site), method = "lm", se = TRUE,
			linewidth = 0.8, alpha = 0.15) +
		# overall trend across all 12 samples
		geom_smooth(method = "lm", se = FALSE, color = "black",
			linetype = "dashed", linewidth = 1) +
		geom_point(aes(color = site, shape = site), size = 3.5) +
		scale_color_manual(values = site_colors,
			labels = c("Cola de Ballena", "Estero Zacatecas")) +
		scale_shape_manual(values = c(16, 17),
			labels = c("Cola de Ballena", "Estero Zacatecas")) +
		labs(x = sed_labels[v], y = "Simpson diversity",
			color = "Site", shape = "Site") +
		theme(legend.position = "bottom",
			legend.title = element_text(size = 9),
			legend.text  = element_text(size = 9))
})

g_scatter <- plot_grid(plotlist = scatter_plots, nrow = 1)
save_plot(
	file.path(output_dir, "SedimentSiteCheck_Simpson_Scatter_February.png"),
	g_scatter, ncol = 3, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


##########################################################################
## 4. Summary printout
##########################################################################

cat("\n=== Wilcoxon tests: sediment variables by site (February) ===\n")
print(wilcox_df, row.names = FALSE)
cat("\nInterpretation guide:\n")
cat("  If p_adj < 0.05 for sed_P, sed_sand or sed_clay: the sites differ\n")
cat("  systematically in those properties — the Feb Spearman correlations\n")
cat("  primarily reflect a site contrast, not a within-site gradient.\n")
cat("  Inspect the scatter plots: if points cluster by colour with a gap\n")
cat("  between clusters, the gradient is site-driven; if colours are\n")
cat("  interleaved along the x-axis, there is a genuine within-site gradient.\n")


##########################################################################
## 5. Full sediment profile: all 9 variables, both sites
##
## Sediment was measured once per transect (n = 6 per site).
## Using unique per-transect values avoids double-counting caused by
## applying the same measurement to both seasons.
##########################################################################

all_sed_vars <- c("sed_pH", "sed_EC", "sed_OM", "sed_orgN",
	"sed_P", "sed_K", "sed_sand", "sed_silt", "sed_clay")

all_sed_labels <- c(
	sed_pH   = "pH",
	sed_EC   = "EC (dS/m)",
	sed_OM   = "OM (%)",
	sed_orgN = "Org-N (%)",
	sed_P    = "P (ppm)",
	sed_K    = "K (ppm)",
	sed_sand = "Sand (%)",
	sed_silt = "Silt (%)",
	sed_clay = "Clay (%)"
)

# Deduplicate: one row per transect (tree-id)
sed_unique <- df_all[!duplicated(df_all[["tree-id"]]),
	c("site", "tree-id", all_sed_vars)]
sed_unique$site <- factor(sed_unique$site, levels = site_order)


## 5a. Wilcoxon tests for all 9 variables -----------------------------------------

wilcox_all <- lapply(all_sed_vars, function(v) {
	cb <- sed_unique[sed_unique$site == "ColaDeBallena",  v]
	ez <- sed_unique[sed_unique$site == "EsteroZacatecas", v]
	wt <- wilcox.test(cb, ez, exact = FALSE)
	data.frame(
		variable  = v,
		CB_median = median(cb, na.rm = TRUE),
		CB_IQR    = IQR(cb,    na.rm = TRUE),
		EZ_median = median(ez, na.rm = TRUE),
		EZ_IQR    = IQR(ez,    na.rm = TRUE),
		W         = wt$statistic,
		p_value   = wt$p.value,
		stringsAsFactors = FALSE
	)
})
wilcox_all_df <- do.call(rbind, wilcox_all)
wilcox_all_df$p_adj <- p.adjust(wilcox_all_df$p_value, method = "BH")
write.csv(wilcox_all_df,
	file.path(output_dir, "SedimentProfile_AllVars_Wilcoxon.csv"),
	row.names = FALSE)

cat("\n=== Wilcoxon tests: all sediment variables by site (unique transects) ===\n")
print(wilcox_all_df, row.names = FALSE)


## 5b. 3×3 boxplot panel — all 9 sediment variables --------------------------------

box_all <- lapply(all_sed_vars, function(v) {
	p_adj   <- wilcox_all_df$p_adj[wilcox_all_df$variable == v]
	y_max   <- max(sed_unique[[v]], na.rm = TRUE)
	y_range <- diff(range(sed_unique[[v]], na.rm = TRUE))

	ggplot(sed_unique, aes(x = site, y = .data[[v]], fill = site)) +
		geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
		geom_jitter(aes(color = site), width = 0.12, size = 2.5, alpha = 0.9) +
		annotate("segment",
			x = 1, xend = 2,
			y = y_max + 0.08 * y_range,
			yend = y_max + 0.08 * y_range) +
		annotate("text",
			x = 1.5, y = y_max + 0.18 * y_range,
			label = sig_label(p_adj), size = 5) +
		scale_fill_manual(values  = site_colors) +
		scale_color_manual(values = site_colors) +
		scale_x_discrete(labels = c("CB", "EZ")) +
		labs(y = all_sed_labels[v], x = NULL,
			title = all_sed_labels[v]) +
		theme(legend.position = "none",
			plot.title  = element_text(size = 10, face = "bold"),
			axis.text.x = element_text(size = 10))
})

g_all_boxes <- plot_grid(plotlist = box_all, nrow = 3, ncol = 3)
save_plot(
	file.path(output_dir, "SedimentProfile_AllVars_Boxplots.png"),
	g_all_boxes, nrow = 3, ncol = 3,
	base_aspect_ratio = 0.9, dpi = 300, base_height = 4)


## 5c. Radar chart — normalized site mean profiles ---------------------------------
##
## Each variable is scaled to [0, 1] using the full observed range across both
## sites, so all axes are comparable. The polygon shows the site mean; individual
## transect points are overlaid on each axis to convey within-site spread.

# Normalize using full range across both sites
mins_all <- apply(sed_unique[, all_sed_vars], 2, min, na.rm = TRUE)
maxs_all <- apply(sed_unique[, all_sed_vars], 2, max, na.rm = TRUE)

norm_fn <- function(x, v) (x - mins_all[v]) / (maxs_all[v] - mins_all[v])

# Site means of normalized values
site_means_norm <- as.data.frame(
	t(sapply(site_order, function(s) {
		rows <- sed_unique[sed_unique$site == s, all_sed_vars]
		sapply(all_sed_vars, function(v) mean(norm_fn(rows[[v]], v), na.rm = TRUE))
	}))
)
site_means_norm$site <- site_order

# Angles: n variables equally spaced, starting at top (pi/2), going clockwise
n_vars  <- length(all_sed_vars)
angles  <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_vars + 1)[1:n_vars]

# Long format for site means (polygon)
radar_mean <- pivot_longer(site_means_norm, cols = -site,
	names_to = "variable", values_to = "value")
radar_mean$angle <- angles[match(radar_mean$variable, all_sed_vars)]
radar_mean$x <- radar_mean$value * cos(radar_mean$angle)
radar_mean$y <- radar_mean$value * sin(radar_mean$angle)

# Close each polygon
radar_closed <- bind_rows(lapply(split(radar_mean, radar_mean$site), function(d) {
	bind_rows(d, d[1, ])
}))

# Individual transect points (normalized)
sed_norm_long <- do.call(rbind, lapply(all_sed_vars, function(v) {
	data.frame(
		site     = sed_unique$site,
		variable = v,
		value    = norm_fn(sed_unique[[v]], v),
		angle    = angles[match(v, all_sed_vars)],
		stringsAsFactors = FALSE
	)
}))
sed_norm_long$x <- sed_norm_long$value * cos(sed_norm_long$angle)
sed_norm_long$y <- sed_norm_long$value * sin(sed_norm_long$angle)
sed_norm_long$site <- factor(sed_norm_long$site, levels = site_order)

# Background grid circles (0.25, 0.5, 0.75, 1.0)
grid_circles <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1.0), function(r) {
	a <- seq(0, 2 * pi, length.out = 200)
	data.frame(x = r * cos(a), y = r * sin(a), r = factor(r))
}))

# Axis spokes and outer variable labels
spoke_df <- data.frame(
	x = cos(angles), y = sin(angles), variable = all_sed_vars)
label_r  <- 1.28
label_df <- data.frame(
	x     = label_r * cos(angles),
	y     = label_r * sin(angles),
	label = unname(all_sed_labels[all_sed_vars]),
	hjust = ifelse(cos(angles) > 0.1, 0, ifelse(cos(angles) < -0.1, 1, 0.5)),
	vjust = ifelse(sin(angles) > 0.1, 0, ifelse(sin(angles) < -0.1, 1, 0.5))
)

radar_plot <- ggplot() +
	# grid circles
	geom_path(data = grid_circles,
		aes(x = x, y = y, group = r),
		color = "grey82", linewidth = 0.35) +
	# axis spokes
	geom_segment(data = spoke_df,
		aes(x = 0, y = 0, xend = x, yend = y),
		color = "grey72", linewidth = 0.35) +
	# site polygons (filled)
	geom_polygon(data = radar_closed,
		aes(x = x, y = y, group = site, fill = site),
		alpha = 0.22) +
	geom_path(data = radar_closed,
		aes(x = x, y = y, group = site, color = site),
		linewidth = 1.1) +
	# individual transect points
	geom_point(data = sed_norm_long,
		aes(x = x, y = y, color = site, shape = site),
		size = 2.2, alpha = 0.85) +
	# grid value labels
	annotate("text",
		x = 0.04, y = c(0.25, 0.5, 0.75, 1.0) + 0.03,
		label = c("0.25", "0.50", "0.75", "1.00"),
		size = 2.6, color = "grey55", hjust = 0) +
	# variable name labels
	geom_text(data = label_df,
		aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
		size = 3.3, lineheight = 0.88) +
	scale_fill_manual(values  = site_colors,
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	scale_color_manual(values = site_colors,
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	scale_shape_manual(values = c(16, 17),
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45)) +
	labs(title = "Sediment fingerprint — Cola de Ballena vs. Estero Zacatecas",
		subtitle = "Polygon = site mean  |  Points = individual transects  |  Axes normalized to [0,1]") +
	theme_void() +
	theme(
		legend.position  = "bottom",
		legend.title     = element_blank(),
		legend.text      = element_text(size = 11),
		plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
		plot.subtitle    = element_text(size = 9,  color = "grey40", hjust = 0.5),
		plot.margin      = margin(10, 10, 10, 10)
	)

save_plot(
	file.path(output_dir, "SedimentProfile_RadarChart.png"),
	radar_plot, base_aspect_ratio = 1.05, dpi = 300, base_height = 8)


##########################################################################
## 6. Combined site fingerprint — intertidal water + sediment
##
## All 15 environmental variables (6 water + 9 sediment) on one radar,
## coloured by site only.
##
## Water variables    → site means averaged across both seasons (n = 12/site)
## Sediment variables → site means from unique transect values   (n =  6/site)
##
## Each variable is independently normalized to [0, 1] so that variables
## with very different units (ppm vs. % vs. °C) are all comparable on
## the same axes.
##########################################################################

# Colours specific to this chart (override earlier site_colors)
fp_site_colors <- c("ColaDeBallena"   = "#e07070",   # light red
	                "EsteroZacatecas" = "#70b870")    # light green
fp_spoke_water  <- "#4393c3"   # blue
fp_spoke_sed    <- "#c9a96e"   # light brown / tan

water_vars_all <- c("water_COND", "water_TDS", "water_ORP",
	"water_pH", "water_salinity", "water_temperature")

combined_vars <- c(water_vars_all, all_sed_vars)

combined_labels <- c(
	water_COND        = "Cond.\n(mS/cm)",
	water_TDS         = "TDS\n(g/L)",
	water_ORP         = "ORP\n(mV)",
	water_pH          = "w-pH",
	water_salinity    = "Salinity\n(PSU)",
	water_temperature = "Temp.\n(°C)",
	sed_pH            = "s-pH",
	sed_EC            = "EC\n(dS/m)",
	sed_OM            = "OM\n(%)",
	sed_orgN          = "Org-N\n(%)",
	sed_P             = "P\n(ppm)",
	sed_K             = "K\n(ppm)",
	sed_sand          = "Sand\n(%)",
	sed_silt          = "Silt\n(%)",
	sed_clay          = "Clay\n(%)"
)

# Site means: water averaged over all 24 samples (both seasons), sediment
# from unique transects — then bind into one row per site
water_site_means <- as.data.frame(
	t(sapply(site_order, function(s) {
		colMeans(df_all[df_all$site == s, water_vars_all], na.rm = TRUE)
	}))
)
sed_site_means <- as.data.frame(
	t(sapply(site_order, function(s) {
		colMeans(sed_unique[sed_unique$site == s, all_sed_vars], na.rm = TRUE)
	}))
)
combined_means <- cbind(water_site_means, sed_site_means)
combined_means$site <- site_order

# Normalize using the full range of INDIVIDUAL data (not site means).
# This ensures site-mean polygon vertices and individual sample points
# both stay within [0, 1], preventing out-of-bounds clipping on the radar.
mins_comb <- c(
	apply(df_all[, water_vars_all], 2, min, na.rm = TRUE),
	apply(sed_unique[, all_sed_vars], 2, min, na.rm = TRUE)
)
maxs_comb <- c(
	apply(df_all[, water_vars_all], 2, max, na.rm = TRUE),
	apply(sed_unique[, all_sed_vars], 2, max, na.rm = TRUE)
)
norm_comb <- function(x, v) (x - mins_comb[v]) / (maxs_comb[v] - mins_comb[v])

# Normalized site means — build as ordered data frame to guarantee spoke order
norm_means_comb <- do.call(rbind, lapply(site_order, function(s) {
	vals <- sapply(combined_vars, function(v)
		norm_comb(combined_means[combined_means$site == s, v], v))
	data.frame(site = s, variable = combined_vars, value = vals,
		row.names = NULL, stringsAsFactors = FALSE)
}))

# Individual sample points: water from df_all, sediment from sed_unique
norm_pts_water <- do.call(rbind, lapply(water_vars_all, function(v) {
	data.frame(site = as.character(df_all$site), variable = v,
		value = norm_comb(df_all[[v]], v), stringsAsFactors = FALSE)
}))
norm_pts_sed <- do.call(rbind, lapply(all_sed_vars, function(v) {
	data.frame(site = as.character(sed_unique$site), variable = v,
		value = norm_comb(sed_unique[[v]], v), stringsAsFactors = FALSE)
}))
norm_pts_all <- rbind(norm_pts_water, norm_pts_sed)
norm_pts_all$site <- factor(norm_pts_all$site, levels = site_order)

# Angles: 15 variables equally spaced, start at top, clockwise
n_comb      <- length(combined_vars)
angles_comb <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_comb + 1)[1:n_comb]
angle_map   <- setNames(angles_comb, combined_vars)

# Coordinates for site-mean polygon — ordered by combined_vars to avoid
# self-intersecting paths (which would make the polygon invisible)
norm_means_comb$angle <- angle_map[norm_means_comb$variable]
norm_means_comb$x     <- norm_means_comb$value * cos(norm_means_comb$angle)
norm_means_comb$y     <- norm_means_comb$value * sin(norm_means_comb$angle)

# Close each polygon by appending first row of each site
radar_closed_comb <- bind_rows(lapply(site_order, function(s) {
	d <- norm_means_comb[norm_means_comb$site == s, ]
	d <- d[order(match(d$variable, combined_vars)), ]
	bind_rows(d, d[1, ])
}))

# Coordinates for individual points
norm_pts_all$angle <- angle_map[norm_pts_all$variable]
norm_pts_all$x     <- norm_pts_all$value * cos(norm_pts_all$angle)
norm_pts_all$y     <- norm_pts_all$value * sin(norm_pts_all$angle)

# Background grid and spokes
grid_comb <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1.0), function(r) {
	a <- seq(0, 2 * pi, length.out = 200)
	data.frame(x = r * cos(a), y = r * sin(a), r = factor(r))
}))
spoke_comb <- data.frame(
	x = cos(angles_comb), y = sin(angles_comb), variable = combined_vars,
	type = c(rep("water", 6), rep("sediment", 9)))

# Variable labels
label_r_comb <- 1.32
label_comb_df <- data.frame(
	x     = label_r_comb * cos(angles_comb),
	y     = label_r_comb * sin(angles_comb),
	label = unname(combined_labels[combined_vars]),
	hjust = ifelse(cos(angles_comb) >  0.1, 0,
		    ifelse(cos(angles_comb) < -0.1, 1, 0.5)),
	vjust = ifelse(sin(angles_comb) >  0.1, 0,
		    ifelse(sin(angles_comb) < -0.1, 1, 0.5))
)

radar_comb_plot <- ggplot() +
	# grid circles
	geom_path(data = grid_comb,
		aes(x = x, y = y, group = r),
		color = "grey82", linewidth = 0.35) +
	# spokes — water in blue, sediment in light brown
	geom_segment(data = spoke_comb[spoke_comb$type == "water", ],
		aes(x = 0, y = 0, xend = x, yend = y),
		color = fp_spoke_water, linewidth = 0.6) +
	geom_segment(data = spoke_comb[spoke_comb$type == "sediment", ],
		aes(x = 0, y = 0, xend = x, yend = y),
		color = fp_spoke_sed, linewidth = 0.6) +
	# site mean polygons — filled shade + border path
	geom_polygon(data = radar_closed_comb,
		aes(x = x, y = y, group = site, fill = site),
		alpha = 0.45, color = NA) +
	geom_path(data = radar_closed_comb,
		aes(x = x, y = y, group = site, color = site),
		linewidth = 1.3) +
	# individual points
	geom_point(data = norm_pts_all,
		aes(x = x, y = y, color = site, shape = site),
		size = 1.8, alpha = 0.75) +
	# grid value labels
	annotate("text",
		x = 0.04, y = c(0.25, 0.5, 0.75, 1.0) + 0.03,
		label = c("0.25", "0.50", "0.75", "1.00"),
		size = 2.4, color = "grey55", hjust = 0) +
	# spoke-type legend as manual annotation
	annotate("segment", x = -1.42, xend = -1.22, y =  1.30, yend =  1.30,
		color = fp_spoke_water, linewidth = 1.4) +
	annotate("text", x = -1.18, y = 1.30, label = "Water",
		size = 3.2, hjust = 0, color = fp_spoke_water) +
	annotate("segment", x = -1.42, xend = -1.22, y =  1.14, yend =  1.14,
		color = fp_spoke_sed, linewidth = 1.4) +
	annotate("text", x = -1.18, y = 1.14, label = "Sediment",
		size = 3.2, hjust = 0, color = fp_spoke_sed) +
	# variable labels
	geom_text(data = label_comb_df,
		aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
		size = 3.0, lineheight = 0.85) +
	scale_fill_manual(values  = fp_site_colors,
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	scale_color_manual(values = fp_site_colors,
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	scale_shape_manual(values = c(16, 17),
		labels = c("Cola de Ballena", "Estero Zacatecas")) +
	coord_fixed(xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55)) +
	labs(
		title    = "Environmental fingerprint — Cola de Ballena vs. Estero Zacatecas",
		subtitle = paste0(
			"Blue spokes: intertidal water (season-averaged)  |  ",
			"Brown spokes: sediment (per transect)\n",
			"Shaded polygon = site mean  |  Points = individual samples/transects  |  ",
			"Axes normalized to [0,1]")
	) +
	theme_void() +
	theme(
		legend.position  = "bottom",
		legend.title     = element_blank(),
		legend.text      = element_text(size = 11),
		plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
		plot.subtitle    = element_text(size = 8.5, color = "grey40", hjust = 0.5,
			lineheight = 1.3),
		plot.margin      = margin(10, 10, 10, 10)
	)

save_plot(
	file.path(output_dir, "CombinedEnvFingerprint_BySite.png"),
	radar_comb_plot, base_aspect_ratio = 1.05, dpi = 300, base_height = 9)
