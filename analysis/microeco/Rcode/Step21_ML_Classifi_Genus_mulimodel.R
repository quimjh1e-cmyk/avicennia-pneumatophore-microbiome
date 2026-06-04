##
## Compare multiple models on classification — Family level
## Thesis adaptation: site as response variable
##
## Original script: Genus level, cal_caretList for simultaneous multi-model comparison.
## Adapted: cal_caretList is not viable for this dataset because:
##   (a) xgbLinear fails with a known xgboost/caret version incompatibility
##       (ALTLIST Set_elt error); cannot be fixed without downgrading xgboost.
##   (b) svmRadial fails with "Variable constant, cannot scale" — zero-variance features
##       arise in training folds when 9-fold CV splits this small dataset too finely.
##   (c) cal_caretList ignores set_trainControl and always uses 9-fold CV internally;
##       with small n this leaves 1-2 test samples per fold, making ROC impossible.
##
## Replacement: individual cal_train calls per model (same pattern as Steps 19-20),
## each with LOOCV and a held-out test set. AUC from cal_ROC(input="pred") is extracted
## for each model and written to a comparison table for the thesis.
##
## Models compared: rf (Random Forest), glmnet (regularized logistic regression), knn (k-NN)
##   - glmnet: handles high-dimensional data natively via L1/L2 regularization
##   - knn: simple distance-based classifier, robust to scaling issues
##   - xgbLinear and svmRadial excluded (see above)
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
tmp_microtable$cal_abund(rel = TRUE)

taxa_level <- "Family"
y_response <- "site"

######################################################
# train and evaluate each model independently with LOOCV + held-out test set

methods_to_compare <- c("rf", "glmnet", "knn")
auc_results <- data.frame(Model = character(), AUC_train = numeric(), AUC_test = numeric(), stringsAsFactors = FALSE)

for (m in methods_to_compare) {
	set.seed(123)
	t_m <- trans_classifier$new(dataset = tmp_microtable, y.response = y_response, x.predictors = taxa_level)
	t_m$cal_split(prop.train = 3/4)

	t_m$set_trainControl(
		method          = "LOOCV",
		classProbs      = TRUE,
		summaryFunction = caret::twoClassSummary
	)
	t_m$cal_train(method = m)
	t_m$cal_predict()

	# AUC on training data
	t_m$cal_ROC(input = "train")
	auc_train <- max(t_m$res_ROC$res_roc$AUC, na.rm = TRUE)

	# AUC on test data
	t_m$cal_ROC(input = "pred")
	auc_test <- max(t_m$res_ROC$res_roc$AUC, na.rm = TRUE)

	auc_results <- rbind(auc_results, data.frame(Model = m, AUC_train = auc_train, AUC_test = auc_test))

	# save ROC curves (test set) per model
	g1 <- t_m$plot_ROC()
	cowplot::save_plot(
		file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_", m, "_ROC_pred.png")),
		g1, base_aspect_ratio = 1.1, dpi = 300, base_height = 6
	)
	# save the object
	save(t_m, file = file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_model_", m, ".RData")))
}

# write model comparison table
write.csv(auc_results, file.path(output_dir, paste0("Classification_", taxa_level, "_", y_response, "_multimodel_AUC_comparison.csv")), row.names = FALSE)
print(auc_results)
