##
## Alpha diversity boxplots — 4 metrics, 4 groups (season x site)
##

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

##########################################################################
# paths
##########################################################################
input_file <- "Output/1.Amplicon/Stage4_AlphaDiversity/AlphaDiv_4metrics.xlsx"
output_dir <- "Output/1.Amplicon/Stage4_AlphaDiversity"

##########################################################################
# load and reshape data
##########################################################################
df <- read_excel(input_file)

# assign group from sample-id prefix
df <- df %>%
  rename(sample_id = `sample-id`) %>%
  mutate(group = case_when(
    grepl("^FCBA", sample_id) ~ "Feb Cola de Ballena",
    grepl("^FEZA", sample_id) ~ "Feb Estero Zacatecas",
    grepl("^SCBA", sample_id) ~ "Sep Cola de Ballena",
    grepl("^SEZA", sample_id) ~ "Sep Estero Zacatecas"
  )) %>%
  mutate(group = factor(group, levels = c(
    "Feb Cola de Ballena",
    "Feb Estero Zacatecas",
    "Sep Cola de Ballena",
    "Sep Estero Zacatecas"
  )))

# pivot to long format
df_long <- df %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson, PD),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson", "PD")))

##########################################################################
# colour palette (season x site)
##########################################################################
group_colors <- c(
  "Feb Cola de Ballena"    = "#4E9AB8",
  "Feb Estero Zacatecas"   = "#B04E1E",
  "Sep Cola de Ballena"    = "#6AAFC2",
  "Sep Estero Zacatecas"   = "#F2B08A"
)

##########################################################################
# custom labeller: show units / range hint in strip
##########################################################################
metric_labels <- c(
  Observed = "Observed ASVs",
  Shannon  = "Shannon index",
  Simpson  = "Simpson index",
  PD = "Faith's PD"
)

##########################################################################
# build plot
##########################################################################
p <- ggplot(df_long, aes(x = group, y = value, fill = group)) +

  # boxplot layer
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85,
               color = "grey30", linewidth = 0.4) +

  # individual data points (jittered)
  geom_jitter(aes(color = group), width = 0.12, size = 2,
              alpha = 0.9, show.legend = FALSE) +

  # one panel per metric, each with its own y-axis
  facet_wrap(~ metric, scales = "free_y", nrow = 2,
             labeller = labeller(metric = metric_labels)) +

  # colours
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = c(
    "Feb Cola de Ballena"    = "#2E7A98",
    "Feb Estero Zacatecas"   = "#B04E1E",
    "Sep Cola de Ballena"    = "#6AAFC2",
    "Sep Estero Zacatecas"   = "#D4845A"
  ), guide = "none") +

  # labels
  labs(x = NULL, y = NULL) +

  # theme
  theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey92", color = "grey60"),
    strip.text        = element_text(face = "bold", size = 11),
    axis.text.x       = element_text(size = 9.5, color = "black", angle = 35, hjust = 1),
    axis.text.y       = element_text(size = 9.5, color = "black"),
    legend.position   = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin       = margin(10, 10, 20, 10)
  )

# save
ggsave(
  filename = file.path(output_dir, "AlphaDiv_boxplots_4groups.png"),
  plot     = p,
  width    = 9, height = 7, dpi = 300
)

message("Saved: AlphaDiv_boxplots_4groups.png")
