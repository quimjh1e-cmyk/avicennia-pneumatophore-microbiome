##
## Regression at Family level
## Thesis adaptation: temperature as continuous response variable
##
## Original script: Genus level, pH as response
## Adapted:
##   Taxonomic level: Genus → Family (36% vs 70% ASV classification rate at genus vs family)
##   Response variable: pH → temperature
##     Rationale: temperature is the key abiotic driver identified in alpha diversity analysis
##     (negative correlation with diversity indices, p < 0.05) and connects ML results
##     directly to the spatiotemporal conclusions of the thesis.
##
## Key methodological changes vs. original protocol:
##   - Taxonomic level: Genus → Family
##   - Response variable: pH → temperature
##   - LOOCV replaces 10-fold CV for single-model training (small n field dataset)
##   - cal_caretList replaced with individual cal_train calls (same fix as Step 21):
##       * cal_caretList ignores set_trainControl and uses 9-fold CV internally → all
##         metrics NA at small n (too few test samples per fold for RMSE computation)
##       * xgbLinear fails with a version incompatibility in the xgboost package
##       * lm fails with p >> n (356 family features, small n → singular matrix)
##     Replacement models: rf, glmnet (L1/L2 regularization handles p >> n), knn
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
# load raw (unnormalized) microtable
input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData"
if(! file.exists(input_path)){
	stop("Please first run the scripts in Stage2 !")
}
load(input_path)
######################################################
# fix a random seed for reproducibility
set.seed(123)

tmp_microtable <- clone(amplicon_16S_microtable)
# calculate relative abundance at each taxonomic level
tmp_microtable$cal_abund(rel = TRUE)


######################################################
# Family level: relative abundance
# Genus level (original) replaced: 36% ASV classification rate at genus vs. 70% at family.
taxa_level <- "Family"

# temperature: continuous response variable
# Negatively correlated with alpha diversity (Spearman r ≈ -0.61, p < 0.05, all samples).
y_response <- "temperature"

######################################################
## Single model: random forest (primary)
######################################################

# create trans_classifier object (handles both classification and regression)
t1 <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)

# split into training and test sets
t1$cal_split(prop.train = 3/4)
write.csv(t1$data_train, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_split_data_train.csv")))
write.csv(t1$data_test,  file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_split_data_test.csv")))

# LOOCV: more appropriate than 10-fold CV for small n
t1$set_trainControl(method = "LOOCV")
t1$cal_train(method = "rf")

# prediction on test set
t1$cal_predict()

# feature importance: IncNodePurity (regression equivalent of MeanDecreaseGini)
t1$cal_feature_imp(scale = TRUE)
write.csv(t1$res_feature_imp, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_feature_importance.csv")))
g1 <- t1$plot_feature_imp(coord_flip = FALSE, colour = "red", fill = "red", width = 0.6) + ylab("IncNodePurity")
cowplot::save_plot(file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_feature_importance.png")), g1,
	base_aspect_ratio = 1.3, dpi = 300, base_height = 6)

save(t1, file = file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_rf_train_predict.RData")))


######################################################
## Multi-model comparison: rf, glmnet, knn
## Individual cal_train calls with the same split and LOOCV — bypasses cal_caretList.
## Test-set RMSE and R² are extracted for each model and written to a comparison table.
######################################################

methods_to_compare <- c("rf", "glmnet", "knn")
regression_results <- data.frame(
	Model      = character(),
	RMSE_LOOCV = numeric(),
	Rsq_LOOCV  = numeric(),
	MAE_LOOCV  = numeric(),
	stringsAsFactors = FALSE
)

for (m in methods_to_compare) {
	set.seed(123)
	t_m <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)
	t_m$cal_split(prop.train = 3/4)
	# savePredictions = "final" stores LOOCV held-out predictions in res_train$pred
	# with standardised "obs" and "pred" columns — avoids data_test column name ambiguity
	t_m$set_trainControl(method = "LOOCV", savePredictions = "final")
	t_m$cal_train(method = m)

	# LOOCV metrics: every training sample held out once — more reliable than the
	# small test set for comparing models at small n
	loocv_preds <- t_m$res_train$pred
	metrics <- caret::postResample(pred = loocv_preds$pred, obs = loocv_preds$obs)

	regression_results <- rbind(regression_results, data.frame(
		Model      = m,
		RMSE_LOOCV = metrics["RMSE"],
		Rsq_LOOCV  = metrics["Rsquared"],
		MAE_LOOCV  = metrics["MAE"]
	))

	save(t_m, file = file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_model_", m, ".RData")))
}

write.csv(regression_results, file.path(output_dir, paste0("Regression_", taxa_level, "_", y_response, "_multimodel_comparison.csv")), row.names = FALSE)
print(regression_results)
