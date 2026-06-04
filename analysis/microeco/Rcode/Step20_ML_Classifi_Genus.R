##
## Classification at Family level
## Thesis adaptation: season (primary) and site (secondary)
##
## Original script: Genus level classification
## Adapted to Family level: only 36% of ASVs were classified to genus level vs. 70% to family
## level, so genus-level feature importance would be dominated by unclassified entries.
## Family level provides named, ecologically interpretable features consistent with the
## differential abundance analysis already conducted at this level.
##
## Original response variables: Cropping (2-group), Fertilization (3-group)
## Adapted:
##   Run 1 — season (February vs September)             [maps to: Cropping]
##   Run 2 — site   (ColaDeBallena vs EsteroZacatecas)  [maps to: Fertilization, now 2-group]
##
## Key methodological changes vs. original protocol:
##   - Taxonomic level: Genus → Family (dataset-specific: 36% vs 70% classification rate)
##   - LOOCV replaces 10-fold CV (more appropriate for small n field datasets)
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
# load raw (unnormalized) microtable; relative abundance computed below
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable)
# calculate relative abundance (total sum scaling) at each taxonomic level
tmp_microtable$cal_abund(rel = TRUE)


######################################################
# Family level: relative abundance
# Genus level (original) replaced: 36% ASV classification rate at genus vs. 70% at family.
taxa_level <- "Family"

######################################################
## Run 1: season — maps to original "Cropping" (2-group single factor)
######################################################
y_response <- "season"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

# split into training and test sets
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

# Boruta feature selection on training data
t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_traindata.csv")))

# LOOCV: more appropriate than 10-fold CV for small n
t1$set_trainControl(method = "LOOCV")
t1$cal_train(method = "rf")

t1$cal_predict()
save(t1, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))

g1 <- t1$plot_confusionMatrix()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_confusionMatrix.png")), g1, base_aspect_ratio = 1.5, dpi = 300, base_height = 6)

# ROC and PR — training data
t1$cal_ROC(input = "train")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_train.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# ROC and PR — test (prediction) data
t1$cal_ROC(input = "pred")
write.csv(t1$res_ROC$res_roc, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.csv")))
write.csv(t1$res_ROC$res_pr,  file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.csv")))
g1 <- t1$plot_ROC()
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_ROC_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)
g1 <- t1$plot_ROC(plot_type = "PR")
cowplot::save_plot(file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_rf_PR_pred.png")), g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6)

# feature importance: MeanDecreaseGini
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
## Run 2: site — maps to original "Fertilization" (now 2-group, not 3-group)
## Site is the dominant structuring factor; feature importance complements Step19 ASV-level results.
######################################################
y_response <- "site"

t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_split_data_test.csv")))

t1$cal_feature_sel(boruta.maxRuns = 300, boruta.pValue = 0.05)
write.csv(t1$data_feature, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_featuresel_boruta_featuredata.csv")))

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
