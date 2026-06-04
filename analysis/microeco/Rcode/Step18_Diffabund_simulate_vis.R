##
## Visualisation of power and FDR from differential abundance benchmarking
## on SparseDOSSA2-simulated data (24 samples, 1000 features, 100 truly differential)
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)

######################################################
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"
output_dir_test <- file.path(output_dir, "2.sim_test")
if(! dir.exists(output_dir_test)){
	stop("Please first run Step 17!")
}


######################################################
# collect results across all methods and simulation runs

all_methods <- list.files(output_dir_test)
res_table   <- data.frame()

for(i in all_methods){
	tmp_list  <- list()
	tmp_path  <- file.path(output_dir_test, i)
	tmp_files <- list.files(tmp_path)
	for(j in tmp_files){
		load(file.path(tmp_path, j))
		tmp_list[[j]] <- c(i, gsub(".RData$", "", j), res[1], res[2])
	}
	tmp       <- do.call(rbind, tmp_list)
	res_table <- rbind(res_table, unname(tmp))
}

colnames(res_table) <- c("method", "param", "power", "fdr")

# parse simulation parameters from file name
# file format: microtable_nsample_24_nfeature_1000_sigfeatures_100_run_k
res_table <- tidyr::separate_wider_delim(res_table, cols = "param",
	names = c("param", "iter"),        delim = "_run_",         cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param",
	names = c("param", "signum"),      delim = "_sigfeatures_", cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param",
	names = c("param", "feature_num"), delim = "_nfeature_",    cols_remove = TRUE)
res_table <- tidyr::separate_wider_delim(res_table, cols = "param",
	names = c("param", "nsample"),     delim = "_nsample_",     cols_remove = TRUE)

# convert character columns to numeric where possible
res_table %<>% dropallfactors(unfac2num = TRUE)

# rename normalization+wilcox methods for clarity
res_table$method[res_table$method == "GMPR"]   <- "GMPR+Wilcox"
res_table$method[res_table$method == "Wrench"] <- "Wrench+Wilcox"


######################################################
# visualise power and FDR across methods

tmp_theme <- theme(
	axis.text.x  = element_text(size = 12),
	axis.text.y  = element_text(size = 10),
	axis.title.y = element_text(size = 14)
)

# power (sensitivity) boxplot
p1 <- res_table %>%
	ggboxplot(x = "method", y = "power", color = "method", palette = "Dark2",
		add = "mean", xlab = "", ylab = "Power (sensitivity)", size = 0.6, width = 0.6)
p1 <- ggpar(p1, legend = "none", x.text.angle = 30) + tmp_theme

# FDR boxplot with reference line at 0.05
p2 <- res_table %>%
	ggboxplot(x = "method", y = "fdr", color = "method", palette = "Dark2",
		add = "mean", xlab = "", ylab = "FDR", size = 0.6, width = 0.6) +
	geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40", linewidth = 0.6)
p2 <- ggpar(p2, legend = "none", x.text.angle = 30) + tmp_theme

# combined figure
g1 <- ggarrange(p1, p2, ncol = 1)
ggexport(g1, filename = file.path(output_dir, "diff_methods_simulation_powerfdr.png"),
	width = 2400, height = 2400, res = 300)

# save summary table for reporting
write.csv(res_table, file.path(output_dir, "diff_methods_simulation_powerfdr.csv"),
	row.names = FALSE)
output_dir <- "./Output/1.Amplicon/Stage6_Diff_abund"

# load one simulation run
load(file.path(output_dir, "1.sim_data",
               "microtable_nsample_24_nfeature_1000_sigfeatures_100_run_1.RData"))

# check group structure
table(tmp_mtobj$sample_table$Group)  # should be 12 vs 12

# verify Feature1-100 are truly differential
# compare mean abundance between groups
grp <- tmp_mtobj$sample_table$Group
f1_grp0 <- mean(as.numeric(tmp_mtobj$otu_table["Feature1", grp == "0"]))
f1_grp1 <- mean(as.numeric(tmp_mtobj$otu_table["Feature1", grp == "1"]))
cat("Feature1 mean Group0:", f1_grp0, "\n")
cat("Feature1 mean Group1:", f1_grp1, "\n")

# compare with a null feature (Feature500, not in spike set)
f500_grp0 <- mean(as.numeric(tmp_mtobj$otu_table["Feature500", grp ==
                                                   "0"]))
f500_grp1 <- mean(as.numeric(tmp_mtobj$otu_table["Feature500", grp ==
                                                   "1"]))
cat("Feature500 mean Group0:", f500_grp0, "\n")
cat("Feature500 mean Group1:", f500_grp1, "\n")