##
## Replot feature importance figures with clean family names (no f__ prefix)
## Loads saved model objects — no retraining needed.
##

library(microeco)
library(magrittr)
library(ggplot2)

output_dir <- "Output/1.Amplicon/Stage7_Machine_learning"

clean_names <- function(t1) {
  rownames(t1$res_feature_imp) <- gsub(".*f__", "", rownames(t1$res_feature_imp))
  t1
}
bigger_y <- theme(axis.text.y = element_text(size = 11))

## --- Classification: Family × site ---
load(file.path(output_dir, "Classification_Family_site_model_rf_train_predict.RData"))

set.seed(123)
t1$cal_feature_imp(scale = TRUE)
clean_names(t1)
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini") + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_site_model_rf_feature_importance.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

set.seed(123)
t1$cal_feature_imp(rf_feature_sig = TRUE, num.rep = 500)
clean_names(t1)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE) + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_site_model_rf_feature_importance_withsig_gini.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE) + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_site_model_rf_feature_importance_withsig_accuracy.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
message("Done: site classification")

## --- Classification: Family × season ---
load(file.path(output_dir, "Classification_Family_season_model_rf_train_predict.RData"))

set.seed(123)
t1$cal_feature_imp(scale = TRUE)
clean_names(t1)
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini") + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_season_model_rf_feature_importance.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

set.seed(123)
t1$cal_feature_imp(rf_feature_sig = TRUE, num.rep = 500)
clean_names(t1)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE) + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_season_model_rf_feature_importance_withsig_gini.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE) + bigger_y
cowplot::save_plot(file.path(output_dir, "Classification_Family_season_model_rf_feature_importance_withsig_accuracy.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
message("Done: season classification")

## --- Regression: Family × temperature ---
load(file.path(output_dir, "Regression_Family_temperature_model_rf_train_predict.RData"))

t1$cal_feature_imp(scale = TRUE)
clean_names(t1)
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("IncNodePurity") + bigger_y
cowplot::save_plot(file.path(output_dir, "Regression_Family_temperature_model_rf_feature_importance.png"), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
message("Done: temperature regression")
