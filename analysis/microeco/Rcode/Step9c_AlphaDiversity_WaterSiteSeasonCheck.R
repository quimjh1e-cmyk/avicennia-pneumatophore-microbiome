##
## Diagnostic: is the water temperature–alpha diversity correlation driven by
##             seasonal separation or a genuine continuous gradient?
##
## Key finding from Step 9 (all samples, water variables):
##   water_temperature shows significant negative Spearman correlation with
##   richness/diversity metrics (Observed, Chao1, ACE, Shannon, Fisher, Pielou)
##   but NOT with evenness metrics (Simpson, InvSimpson, se.chao1, Coverage).
##
## Unlike sediment (constant across seasons), water variables vary with BOTH
## season and site. Comparisons are therefore done along both dimensions.
##
## Statistical test: Wilcoxon rank-sum, non-parametric, n = 12 per group.
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
meta     <- tmp_microtable$sample_table
alphadiv <- tmp_microtable$alpha_diversity
df_all   <- cbind(meta, alphadiv[rownames(meta), ])

# Factor ordering
site_order   <- c("ColaDeBallena", "EsteroZacatecas")
season_order <- c("February", "September")
df_all$site   <- factor(df_all$site,   levels = site_order)
df_all$season <- factor(df_all$season, levels = season_order)

# Colour palettes
site_colors   <- c("ColaDeBallena"   = "#1f78b4", "EsteroZacatecas" = "#e66101")
season_colors <- c("February" = "#2166ac",         "September"      = "#d73027")
season_shapes <- c("February" = 16,                "September"      = 17)

# 4 season×site group colours for radar chart
group_colors <- c(
	"February_ColaDeBallena"    = "#2166ac",
	"February_EsteroZacatecas"  = "#74add1",
	"September_ColaDeBallena"   = "#d73027",
	"September_EsteroZacatecas" = "#f46d43"
)
group_shapes <- c(
	"February_ColaDeBallena"    = 16,
	"February_EsteroZacatecas"  = 17,
	"September_ColaDeBallena"   = 15,
	"September_EsteroZacatecas" = 18
)

# Water variables and display labels (prefix stripped for plots)
water_vars <- c("water_COND", "water_TDS", "water_ORP",
	"water_pH", "water_salinity", "water_temperature")
water_labels <- c(
	water_COND        = "Conductivity\n(mS/cm)",
	water_TDS         = "TDS (g/L)",
	water_ORP         = "ORP (mV)",
	water_pH          = "pH",
	water_salinity    = "Salinity (PSU)",
	water_temperature = "Temperature (°C)"
)

# Significance label helper
sig_label <- function(p) {
	if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else "ns"
}

# Wilcoxon helper: group A vs group B within df
run_wilcox <- function(df, vars, group_col, level_a, level_b) {
	do.call(rbind, lapply(vars, function(v) {
		a <- df[df[[group_col]] == level_a, v]
		b <- df[df[[group_col]] == level_b, v]
		wt <- wilcox.test(a, b, exact = FALSE)
		data.frame(
			variable    = v,
			A_median    = median(a, na.rm = TRUE),
			A_IQR       = IQR(a,    na.rm = TRUE),
			B_median    = median(b, na.rm = TRUE),
			B_IQR       = IQR(b,    na.rm = TRUE),
			W           = wt$statistic,
			p_value     = wt$p.value,
			stringsAsFactors = FALSE
		)
	}))
}


##########################################################################
## 1. Wilcoxon tests — water variables by season (February vs September)
##########################################################################

wilcox_season <- run_wilcox(df_all, water_vars, "season", "February", "September")
colnames(wilcox_season)[2:5] <- c("Feb_median", "Feb_IQR", "Sep_median", "Sep_IQR")
wilcox_season$p_adj <- p.adjust(wilcox_season$p_value, method = "BH")
write.csv(wilcox_season,
	file.path(output_dir, "WaterSiteCheck_Wilcoxon_bySeason.csv"),
	row.names = FALSE)


##########################################################################
## 2. Wilcoxon tests — water variables by site (CB vs EZ)
##########################################################################

wilcox_site <- run_wilcox(df_all, water_vars, "site", "ColaDeBallena", "EsteroZacatecas")
colnames(wilcox_site)[2:5] <- c("CB_median", "CB_IQR", "EZ_median", "EZ_IQR")
wilcox_site$p_adj <- p.adjust(wilcox_site$p_value, method = "BH")
write.csv(wilcox_site,
	file.path(output_dir, "WaterSiteCheck_Wilcoxon_bySite.csv"),
	row.names = FALSE)


##########################################################################
## 3. Scatter plots — Shannon vs. water_temperature (the significant variable)
##
## Two versions:
##   3a. Coloured by season — diagnostic: is the correlation a step-function
##       between seasons or a continuous gradient spanning both?
##   3b. Coloured by site   — diagnostic: does one site drive the signal?
##
## Each includes per-group trend lines (coloured) and an overall dashed line.
##########################################################################

make_scatter <- function(color_var, color_pal, shape_pal, legend_labels,
	filename) {
	ggplot(df_all, aes(x = water_temperature, y = Shannon)) +
		geom_smooth(aes(color = .data[[color_var]]),
			method = "lm", se = TRUE, linewidth = 0.85, alpha = 0.12) +
		geom_smooth(method = "lm", se = FALSE,
			color = "black", linetype = "dashed", linewidth = 1) +
		geom_point(aes(color = .data[[color_var]],
			shape = .data[[color_var]]), size = 3.5) +
		scale_color_manual(values = color_pal, labels = legend_labels) +
		scale_shape_manual(values = shape_pal,  labels = legend_labels) +
		labs(x = "Water temperature (°C)", y = "Shannon diversity",
			color = NULL, shape = NULL) +
		theme(legend.position = "bottom",
			legend.text = element_text(size = 10))
}

g3a <- make_scatter(
	color_var     = "season",
	color_pal     = season_colors,
	shape_pal     = season_shapes,
	legend_labels = c("February", "September"),
	filename      = NULL
)
g3b <- make_scatter(
	color_var     = "site",
	color_pal     = site_colors,
	shape_pal     = c(16, 17),
	legend_labels = c("Cola de Ballena", "Estero Zacatecas"),
	filename      = NULL
)

g_scatter_temp <- plot_grid(g3a, g3b, nrow = 1,
	labels = c(" ", "By site"), label_size = 11,
	label_fontface = "plain", hjust = 0)
save_plot(
	file.path(output_dir, "WaterSiteCheck_Shannon_vs_Temperature.png"),
	g_scatter_temp, ncol = 2, base_aspect_ratio = 1.1, dpi = 300, base_height = 5)


##########################################################################
## 4. Boxplots — all 6 water variables by season (2×3 grid)
##########################################################################

make_boxplot_group <- function(df, group_col, group_levels, color_pal,
	axis_labels, wilcox_df) {

	lapply(water_vars, function(v) {
		p_adj   <- wilcox_df$p_adj[wilcox_df$variable == v]
		y_max   <- max(df[[v]], na.rm = TRUE)
		y_range <- diff(range(df[[v]], na.rm = TRUE))

		ggplot(df, aes(x = .data[[group_col]], y = .data[[v]],
			fill = .data[[group_col]])) +
			geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
			geom_jitter(aes(color = .data[[group_col]]),
				width = 0.12, size = 2.5, alpha = 0.9) +
			annotate("segment",
				x = 1, xend = 2,
				y = y_max + 0.08 * y_range,
				yend = y_max + 0.08 * y_range) +
			annotate("text",
				x = 1.5, y = y_max + 0.18 * y_range,
				label = sig_label(p_adj), size = 5) +
			scale_fill_manual(values  = color_pal) +
			scale_color_manual(values = color_pal) +
			scale_x_discrete(labels = axis_labels) +
			labs(y = water_labels[v], x = NULL, title = water_labels[v]) +
			theme(legend.position = "none",
				plot.title  = element_text(size = 10, face = "bold"),
				axis.text.x = element_text(size = 10))
	})
}

# By season
box_season <- make_boxplot_group(
	df          = df_all,
	group_col   = "season",
	group_levels= season_order,
	color_pal   = season_colors,
	axis_labels = c("February", "September"),
	wilcox_df   = wilcox_season
)
g_box_season <- plot_grid(plotlist = box_season, nrow = 2, ncol = 3)
save_plot(
	file.path(output_dir, "WaterProfile_Boxplots_bySeason.png"),
	g_box_season, nrow = 2, ncol = 3,
	base_aspect_ratio = 0.9, dpi = 300, base_height = 4)

# By site
box_site <- make_boxplot_group(
	df          = df_all,
	group_col   = "site",
	group_levels= site_order,
	color_pal   = site_colors,
	axis_labels = c("CB", "EZ"),
	wilcox_df   = wilcox_site
)
g_box_site <- plot_grid(plotlist = box_site, nrow = 2, ncol = 3)
save_plot(
	file.path(output_dir, "WaterProfile_Boxplots_bySite.png"),
	g_box_site, nrow = 2, ncol = 3,
	base_aspect_ratio = 0.9, dpi = 300, base_height = 4)


##########################################################################
## 5. Radar chart — water fingerprint by season × site (4 groups)
##
## Water conditions vary with both season and site, so showing all 4
## combinations reveals whether temperature separates seasons, sites, or both.
## Each variable is normalized to [0,1] using the full observed range.
##########################################################################

df_all$group <- paste(df_all$season, df_all$site, sep = "_")
group_order  <- names(group_colors)

mins_w <- apply(df_all[, water_vars], 2, min, na.rm = TRUE)
maxs_w <- apply(df_all[, water_vars], 2, max, na.rm = TRUE)
norm_w <- function(x, v) (x - mins_w[v]) / (maxs_w[v] - mins_w[v])

# Group means of normalized values
group_means_norm <- as.data.frame(
	t(sapply(group_order, function(g) {
		rows <- df_all[df_all$group == g, water_vars]
		sapply(water_vars, function(v) mean(norm_w(rows[[v]], v), na.rm = TRUE))
	}))
)
group_means_norm$group <- group_order

# Angles: 6 variables, start at top, clockwise
n_vars_w <- length(water_vars)
angles_w <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_vars_w + 1)[1:n_vars_w]

# Long format for group means (polygon)
radar_mean_w <- pivot_longer(group_means_norm, cols = -group,
	names_to = "variable", values_to = "value")
radar_mean_w$angle <- angles_w[match(radar_mean_w$variable, water_vars)]
radar_mean_w$x <- radar_mean_w$value * cos(radar_mean_w$angle)
radar_mean_w$y <- radar_mean_w$value * sin(radar_mean_w$angle)
radar_closed_w <- bind_rows(lapply(split(radar_mean_w, radar_mean_w$group), function(d) {
	bind_rows(d, d[1, ])
}))

# Individual sample points (normalized)
samp_norm_long <- do.call(rbind, lapply(water_vars, function(v) {
	data.frame(
		group    = df_all$group,
		variable = v,
		value    = norm_w(df_all[[v]], v),
		angle    = angles_w[match(v, water_vars)],
		stringsAsFactors = FALSE
	)
}))
samp_norm_long$x <- samp_norm_long$value * cos(samp_norm_long$angle)
samp_norm_long$y <- samp_norm_long$value * sin(samp_norm_long$angle)
samp_norm_long$group <- factor(samp_norm_long$group, levels = group_order)

# Background grid and spokes
grid_circles_w <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1.0), function(r) {
	a <- seq(0, 2 * pi, length.out = 200)
	data.frame(x = r * cos(a), y = r * sin(a), r = factor(r))
}))
spoke_w <- data.frame(x = cos(angles_w), y = sin(angles_w), variable = water_vars)

# Strip water_ prefix for radar labels
water_labels_short <- c(
	water_COND        = "Cond.\n(mS/cm)",
	water_TDS         = "TDS\n(g/L)",
	water_ORP         = "ORP\n(mV)",
	water_pH          = "pH",
	water_salinity    = "Salinity\n(PSU)",
	water_temperature = "Temp.\n(°C)"
)
label_r_w <- 1.30
label_w_df <- data.frame(
	x     = label_r_w * cos(angles_w),
	y     = label_r_w * sin(angles_w),
	label = unname(water_labels_short[water_vars]),
	hjust = ifelse(cos(angles_w) > 0.1, 0, ifelse(cos(angles_w) < -0.1, 1, 0.5)),
	vjust = ifelse(sin(angles_w) > 0.1, 0, ifelse(sin(angles_w) < -0.1, 1, 0.5))
)

legend_labels_w <- c(
	"February_ColaDeBallena"    = "February — Cola de Ballena",
	"February_EsteroZacatecas"  = "February — Estero Zacatecas",
	"September_ColaDeBallena"   = "September — Cola de Ballena",
	"September_EsteroZacatecas" = "September — Estero Zacatecas"
)

radar_plot_w <- ggplot() +
	geom_path(data = grid_circles_w,
		aes(x = x, y = y, group = r),
		color = "grey82", linewidth = 0.35) +
	geom_segment(data = spoke_w,
		aes(x = 0, y = 0, xend = x, yend = y),
		color = "grey72", linewidth = 0.35) +
	geom_polygon(data = radar_closed_w,
		aes(x = x, y = y, group = group, fill = group),
		alpha = 0.18) +
	geom_path(data = radar_closed_w,
		aes(x = x, y = y, group = group, color = group),
		linewidth = 1.0) +
	geom_point(data = samp_norm_long,
		aes(x = x, y = y, color = group, shape = group),
		size = 2.0, alpha = 0.80) +
	annotate("text",
		x = 0.04, y = c(0.25, 0.5, 0.75, 1.0) + 0.03,
		label = c("0.25", "0.50", "0.75", "1.00"),
		size = 2.5, color = "grey55", hjust = 0) +
	geom_text(data = label_w_df,
		aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
		size = 3.3, lineheight = 0.88) +
	scale_fill_manual(values  = group_colors, labels = legend_labels_w) +
	scale_color_manual(values = group_colors, labels = legend_labels_w) +
	scale_shape_manual(values = group_shapes,  labels = legend_labels_w) +
	coord_fixed(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45)) +
	labs(
		title    = "Intertidal water fingerprint — season × site",
		subtitle = "Polygon = group mean  |  Points = individual samples  |  Axes normalized to [0,1]"
	) +
	theme_void() +
	theme(
		legend.position  = "bottom",
		legend.title     = element_blank(),
		legend.text      = element_text(size = 10),
		legend.key.size  = unit(0.5, "cm"),
		plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
		plot.subtitle    = element_text(size = 9,  color = "grey40", hjust = 0.5),
		plot.margin      = margin(10, 10, 10, 10)
	) +
	guides(color = guide_legend(nrow = 2))

save_plot(
	file.path(output_dir, "WaterProfile_RadarChart_SeasonSite.png"),
	radar_plot_w, base_aspect_ratio = 1.05, dpi = 300, base_height = 8)


##########################################################################
## 6. Summary printout
##########################################################################

cat("\n=== Wilcoxon tests: water variables by season ===\n")
print(wilcox_season, row.names = FALSE)
cat("\n=== Wilcoxon tests: water variables by site ===\n")
print(wilcox_site, row.names = FALSE)
cat("\nInterpretation guide:\n")
cat("  Scatter plot 3a: if the two season clouds are separated along the\n")
cat("  x-axis (temperature) with little overlap, the negative correlation\n")
cat("  between temperature and Shannon is a seasonal step-effect, not a\n")
cat("  continuous dose-response. If the clouds overlap and the within-season\n")
cat("  lines also slope negatively, the relationship is a genuine gradient.\n")
cat("\n")
cat("  Radar chart: if the September polygons (red tones) and February\n")
cat("  polygons (blue tones) are clearly separated but CB/EZ shapes within\n")
cat("  each season overlap, temperature is a seasonal driver with no\n")
cat("  site-specific effect. If CB and EZ polygons diverge within a season,\n")
cat("  site-level differences in water chemistry are also relevant.\n")
