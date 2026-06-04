##
## Replot site feature importance: true top 10 by MeanDecreaseAccuracy
## The default plot_feature_imp selects families by an internal criterion
## that does not match the Accuracy ranking. This script extracts the
## rfPermute results and plots the correct top 10.
##

library(microeco)
library(magrittr)
library(ggplot2)

output_dir <- "Output/1.Amplicon/Stage7_Machine_learning"

load(file.path(output_dir, "Classification_Family_site_model_rf_train_predict.RData"))

set.seed(123)
t1$cal_feature_imp(rf_feature_sig = TRUE, num.rep = 500)
imp <- t1$res_feature_imp

imp$Family <- gsub(".*f__", "", rownames(imp))

imp_sorted <- imp[order(-imp$MeanDecreaseAccuracy), ]
top10 <- head(imp_sorted, 10)

sig_labels <- cut(top10$MeanDecreaseAccuracy.pval,
                  breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("**", "*", "ns"))

top10$Significance <- sig_labels
top10$Family <- factor(top10$Family, levels = rev(top10$Family))

pal <- c("**" = "#1b9e77", "*" = "#d95f02", "ns" = "#7570b3")

g1 <- ggplot(top10, aes(x = Family, y = MeanDecreaseAccuracy, fill = Significance)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = as.character(Significance)),
            vjust = -0.5, size = 4, colour = "grey30") +
  scale_fill_manual(values = pal, name = "Significance") +
  labs(x = NULL, y = "MeanDecreaseAccuracy") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 11, angle = 40, hjust = 1),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

cowplot::save_plot(
  file.path(output_dir, "Classification_Family_site_model_rf_top10_accuracy.png"),
  g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

## Families with MeanDecreaseAccuracy >= 4
top_acc <- imp_sorted[imp_sorted$MeanDecreaseAccuracy >= 4, ]
top_acc$sig_label <- cut(top_acc$MeanDecreaseAccuracy.pval,
                         breaks = c(-Inf, 0.01, 0.05, Inf),
                         labels = c("**", "*", "ns"))
top_acc$Family <- factor(top_acc$Family, levels = rev(top_acc$Family))

g2 <- ggplot(top_acc, aes(x = Family, y = MeanDecreaseAccuracy, fill = sig_label)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = as.character(sig_label)),
            vjust = -0.5, size = 4, colour = "grey30") +
  scale_fill_manual(values = pal, name = "Significance") +
  labs(x = NULL, y = "MeanDecreaseAccuracy") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 11, angle = 40, hjust = 1),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

cowplot::save_plot(
  file.path(output_dir, "Classification_Family_site_model_rf_acc_ge4.png"),
  g2, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

message("Saved: Classification_Family_site_model_rf_top10_accuracy.png")
message("\nFull ranking by MeanDecreaseAccuracy:")
n_sig <- sum(imp_sorted$MeanDecreaseAccuracy.pval < 0.05)
for (i in seq_len(nrow(imp_sorted))) {
  message(sprintf("  %2d. %-30s  Acc=%.4f  p=%.4f",
    i, imp_sorted$Family[i], imp_sorted$MeanDecreaseAccuracy[i],
    imp_sorted$MeanDecreaseAccuracy.pval[i]))
}
message(sprintf("\nSignificant families (p < 0.05): %d / %d", n_sig, nrow(imp_sorted)))
