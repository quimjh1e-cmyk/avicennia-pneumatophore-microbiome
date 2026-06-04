##
## Classification at ASV level
## Thesis adaptation: site (primary) and season (secondary)
##
## Original response variable: Cropping
## Adapted:
##   Run 1 — site   (ColaDeBallena vs EsteroZacatecas)  [primary: strongest PERMANOVA predictor]
##   Run 2 — season (February vs September)             [secondary: weaker signal, useful contrast]
##
## Key methodological changes vs. original protocol:
##   - Boruta feature selection enabled (was commented out)
##   - LOOCV replaces 10-fold CV (more appropriate for small n field datasets)
##   - y_response variable used throughout for consistent file naming
##   - No rhizosphere/compartment subset: full microtable used
##


######################################################
# load packages
library(microeco)
library(magrittr)
######################################################
# create an output directory if it does not exist
output_dir <- "Output/1.Amplicon/Stage7_Machine_learning"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load CLR-normalized microtable
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_CLR.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable_CLR)

# Filter low-abundance ASVs using raw (unnormalized) data, then apply filter to CLR data
load("Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData")
tmp <- clone(amplicon_16S_microtable)
tmp$filter_taxa(rel_abund = 0.001)
tmp_microtable$tax_table %<>% .[rownames(.) %in% rownames(tmp$tax_table), ]
tmp_microtable$tidy_dataset()

# Add ASV identifiers to taxonomy; regenerate taxa_abund using CLR values (rel = FALSE)
tmp_microtable$add_rownames2taxonomy(use_name = "ASV")
tmp_microtable$cal_abund(rel = FALSE)

taxa_level <- "ASV"


######################################################
## Run 1: site — primary classification target
## Site is the strongest predictor of community composition (PERMANOVA).
## Feature importance here triangulates with differential abundance results.
######################################################
y_response <- "site"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

# split data into training and testing sets
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

# Boruta feature selection on training data (Boruta method, 300 max runs)
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_traindata.csv")))

# LOOCV: leave-one-out cross-validation, more appropriate than 10-fold CV for small n
t1$set_trainControl(method = "LOOCV")
t1$cal_train(method = "rf")

# prediction on test set
t1$cal_predict()
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# ROC and PR curves — training data
t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# ROC and PR curves — test (prediction) data
t1$cal_ROC(input = "pred")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# feature importance: MeanDecreaseGini (varImp via caret)
t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

# feature importance with permutation-based significance (rfPermute)
t1$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1,
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1,
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)


######################################################
## Run 2: season — secondary classification target
## Expected lower performance than site (season is a weaker structuring factor).
## Lower AUC here vs. site reinforces the PERMANOVA finding that site > season.
######################################################
y_response <- "season"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_traindata.csv")))

t1$set_trainControl(method = "LOOCV")
t1$cal_train(method = "rf")

t1$cal_predict()
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

t1$cal_ROC(input = "pred")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("MeanDecreaseGini")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1, base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

t1$cal_feature_imp(rf_feature_sig = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig.csv")))
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseGini", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_gini.png")), g1,
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
g1 <- t1$plot_feature_imp(rf_sig_show = "MeanDecreaseAccuracy", show_sig_group = TRUE, coord_flip = FALSE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_feature_importance_withsig_accuracy.png")), g1,
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)
