##
## Simulate microbial communities for benchmarking differential abundance methods
##
## SparseDOSSA2 fits parameters from the real dataset and generates synthetic
## communities with known ground truth (controlled significant feature ratio
## and effect size), allowing evaluation of method performance.
##
## Template: actual ASV-level data (family-filtered, 1346 features, 24 samples)
## Simulation design: binary group (2 x 12 samples), mirroring the season and
## site comparisons in Steps 12-13.
##


######################################################
# load packages
library(SparseDOSSA2)
library(microeco)
library(magrittr)

######################################################
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}

# load preprocessed microtable (family-filtered, consistent with Steps 12-14)
load(file.path(output_dir, "tmp_microtable_preproc.RData"))

# tidy dataset to ensure clean state before fitting
tmp_microtable$tidy_dataset()

######################################################
# fit SparseDOSSA2 parameters from the real ASV data
# this captures the sparsity, compositionality, and
# covariance structure of the actual pneumatophore communities
fitted_param <- fit_SparseDOSSA2(data = as.matrix(tmp_microtable$otu_table))
save(fitted_param,
	file = file.path(output_dir, "microtable_SparseDOSSA2_fitted_param.RData"),
	compress = TRUE)

######################################################
# simulation parameters
# n_sample    : matches actual study design (24 samples, 12 per group)
# n_feature   : close to actual filtered ASV count (~1346)
# sig ratio   : 10 % of features are truly differentially abundant
# effect_size : log fold change of 2 (moderate-strong signal)
# n_sim       : 10 replicate simulations per parameter set

dir_simdata <- file.path(output_dir, "1.sim_data")
dir.create(dir_simdata, showWarnings = FALSE)

n_sample          <- 24
n_feature         <- 1000
n_sigfeature_ratio <- 0.1
n_sigfeature      <- n_sigfeature_ratio * n_feature
n_sim             <- 10

load(file.path(output_dir, "microtable_SparseDOSSA2_fitted_param.RData"))

for(i in n_sample){
	for(j in n_sigfeature){
		group_sample_num <- i / 2
		# binary group: mirrors 2-group design (season or site)
		Group <- c(rep(1, group_sample_num), rep(0, group_sample_num))
		metadata_matrix <- as.matrix(data.frame(Group))
		# spike metadata: first j features are truly differential with effect size 2
		spike_metadata <- data.frame(
			metadata_datum    = 1,
			feature_spiked    = paste0("Feature", 1:j),
			associated_property = "abundance",
			effect_size       = 2
		)

		for(k in seq_len(n_sim)){
			set.seed(k)
			sim1 <- SparseDOSSA2(
				template        = fitted_param,
				n_sample        = i,
				n_feature       = n_feature,
				new_features    = TRUE,
				spike_metadata  = spike_metadata,
				metadata_matrix = metadata_matrix
			)

			tmp_sample <- sim1$spike_metadata$metadata_matrix %>% as.data.frame
			tmp_abund  <- sim1$simulated_data %>% as.data.frame
			rownames(tmp_sample) <- colnames(tmp_abund)
			tmp_tax <- data.frame(ASV = rownames(tmp_abund))
			rownames(tmp_tax) <- tmp_tax[, 1]

			# build microtable object
			tmp_mtobj <- microtable$new(
				sample_table = tmp_sample,
				otu_table    = tmp_abund,
				tax_table    = tmp_tax
			)
			tmp_mtobj$tidy_dataset()
			tmp_mtobj$raw_sim <- sim1

			save(tmp_mtobj,
				file = file.path(dir_simdata,
					paste0("microtable_nsample_", i,
					       "_nfeature_", n_feature,
					       "_sigfeatures_", j,
					       "_run_", k, ".RData")),
				compress = TRUE)
		}
	}
}
