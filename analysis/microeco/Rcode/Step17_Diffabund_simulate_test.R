##
## Differential abundance test on SparseDOSSA2-simulated data
## Evaluates power (sensitivity) and FDR for each method using known ground truth
##
## Methods tested: metagenomeSeq, linda, ALDEx2_kw, ancombc2, DESeq2, edgeR,
##                 GMPR + wilcox, Wrench + wilcox
## Simulated design: binary group (2 x 12), 1000 features, 100 truly differential
##


######################################################
# load packages
library(microeco)
library(magrittr)

######################################################
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	stop("Please first run Step 16!")
}

dir_simdata     <- file.path(output_dir, "1.sim_data")
output_dir_test <- file.path(output_dir, "2.sim_test")
dir.create(output_dir_test, recursive = TRUE, showWarnings = FALSE)


######################################################
# helper: calculate power (sensitivity) and FDR
# H    : ground-truth vector (1 = truly differential, 0 = null)
# test : method result vector (1 = called significant, 0 = not)
powerfdr <- function(H, test){
	tp     <- sum(H != 0 & test != 0)
	fp     <- sum(H == 0 & test != 0)
	fn     <- sum(H != 0 & test == 0)
	power1 <- tp / (tp + fn)
	fdr1   <- fp / (tp + fp)
	c(power1, fdr1)
}


######################################################
# run differential tests on all simulated datasets

all_input  <- list.files(dir_simdata)
group      <- "Group"
taxa_level <- "ASV"

for(method in c("metagenomeSeq", "linda", "ALDEx2_kw", "ancombc2", "DESeq2", "edgeR", "GMPR", "Wrench")){
	dir_difftest_method <- file.path(output_dir_test, method)
	dir.create(dir_difftest_method, showWarnings = FALSE)

	for(i in all_input){
		load(file.path(dir_simdata, i))
		# ensure Group is treated as a categorical variable
		tmp_mtobj$sample_table$Group %<>% as.character

		if(method %in% c("GMPR", "Wrench")){
			tmp      <- trans_norm$new(tmp_mtobj)
			tmp_mtobj <- tmp$norm(method = method, condition = group)
			tmp_mtobj$cal_abund(rel = FALSE)
			tmp_diff <- trans_diff$new(dataset = tmp_mtobj, method = "wilcox",
				group = group, taxa_level = taxa_level)
		} else {
			tmp_diff <- trans_diff$new(dataset = tmp_mtobj, method = method,
				group = group, taxa_level = taxa_level)
		}

		# extract indices of features called significant
		tmp  <- tmp_diff$res_diff %>% .[grepl("*", .$Significance, fixed = TRUE), ] %>% .[, "Taxa"]
		rej  <- gsub("Feature", "", tmp) %>% as.numeric

		# parse simulation parameters from file name
		# file format: microtable_nsample_24_nfeature_1000_sigfeatures_100_run_k.RData
		n_feature <- gsub(".*nfeature_(\\d+)_.*",    "\\1", i) %>% as.numeric
		j         <- gsub(".*sigfeatures_(\\d+)_.*", "\\1", i) %>% as.numeric

		# build ground-truth (H) and test result vectors
		H    <- test <- rep(0, n_feature)
		H[1:j]   <- 1
		test[rej] <- 1

		# calculate and save power + FDR for this simulation run
		res <- powerfdr(H, test)
		save(res, file = file.path(dir_difftest_method, i), compress = TRUE)
	}
}
