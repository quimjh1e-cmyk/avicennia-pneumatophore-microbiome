
### TASK 1 - IMPORT PICRUSt2 OUTPUT INTO MICROECO


#   Step 1 — Load dependencies and check the source data

# Load the rarefied microtable (for sample metadata)
library(microeco)
library(magrittr)

input_path <- "Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
if (!file.exists(input_path)) stop("File not found — check working
  directory with getwd()")
load(input_path)

sample_meta <- amplicon_16S_microtable_rarefy$sample_table

cat("=== Sample metadata ===\n")
cat("Samples:", nrow(sample_meta), "\n")
cat("Columns:", paste(colnames(sample_meta), collapse = ", "), "\n\n")
cat("Season levels:", paste(unique(sample_meta$season), collapse = ", "),
    "\n")
cat("Site levels:  ", paste(unique(sample_meta$site),   collapse = ", "),
    "\n")
print(sample_meta[, c("season", "site")])

# Step 2 — Import the three PICRUSt2 tables

base_path <- "../../picrust2/output"

# --- MetaCyc pathways ---
pw_raw  <- read.delim(file.path(base_path, "pathways_described.tsv"),
                      row.names = 1, check.names = FALSE)
pw_desc <- pw_raw[["description"]]
names(pw_desc) <- rownames(pw_raw)
pw_abun <- pw_raw[, colnames(pw_raw) != "description"]

# --- KO ---
ko_raw  <- read.delim(file.path(base_path, "KO_described.tsv"),
                      row.names = 1, check.names = FALSE)
ko_desc <- ko_raw[["description"]]
names(ko_desc) <- rownames(ko_raw)
ko_abun <- ko_raw[, colnames(ko_raw) != "description"]

# --- EC ---
ec_raw  <- read.delim(file.path(base_path, "EC_described.tsv"),
                      row.names = 1, check.names = FALSE)
ec_desc <- ec_raw[["description"]]
names(ec_desc) <- rownames(ec_raw)
ec_abun <- ec_raw[, colnames(ec_raw) != "description"]

cat("Pathways:", nrow(pw_abun), "x", ncol(pw_abun), "\n")
cat("KOs:     ", nrow(ko_abun), "x", ncol(ko_abun), "\n")
cat("ECs:     ", nrow(ec_abun), "x", ncol(ec_abun), "\n")

common <- intersect(colnames(pw_abun), rownames(sample_meta))
cat("\nSamples in common:", length(common), "\n")
cat("Missing from PICRUSt2:", setdiff(rownames(sample_meta),
                                      colnames(pw_abun)), "\n")
cat("Missing from metadata:", setdiff(colnames(pw_abun),
                                      rownames(sample_meta)), "\n")

## Step 3 — Build microtable objects

pw_abun <- pw_abun[, common]
ko_abun <- ko_abun[, common]
ec_abun <- ec_abun[, common]
meta_sub <- sample_meta[common, ]

func_pw <- microtable$new(otu_table = pw_abun, sample_table = meta_sub)
func_pw$tidy_dataset()

func_ko <- microtable$new(otu_table = ko_abun, sample_table = meta_sub)
func_ko$tidy_dataset()

func_ec <- microtable$new(otu_table = ec_abun, sample_table = meta_sub)
func_ec$tidy_dataset()

cat("func_pw:", nrow(func_pw$otu_table), "pathways x",
    nrow(func_pw$sample_table), "samples\n")
cat("func_ko:", nrow(func_ko$otu_table), "KOs x",
    nrow(func_ko$sample_table), "samples\n")
cat("func_ec:", nrow(func_ec$otu_table), "ECs x",
    nrow(func_ec$sample_table), "samples\n")


  ## Step 4 — Save checkpoint

output_dir <- "Output/4.PICRUSt2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save(func_pw, pw_desc,
     func_ko, ko_desc,
     func_ec, ec_desc,
     file = file.path(output_dir, "picrust2_microtables.RData"))

cat("Saved to:", file.path(output_dir, "picrust2_microtables.RData"),
    "\n")



## Step 5 — Differential Abundance (DA) with 8 Methods for Season, Site, and Season_Site ##

## Step 5 — Differential Abundance (DA) with 8 Methods for Season, Site, and Season_Site ##

# Define output directory
output_dir <- "Output/4.PICRUSt2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create the interaction variable 'season_site'
func_pw$sample_table$season_site <- paste(
  func_pw$sample_table$season, 
  func_pw$sample_table$site, 
  sep = "_"
)

cat("\n=== Applying Integer Fix for Count-Based Methods ===\n")
# Integer copy for count-based methods & normalizations
func_pw_int <- func_pw$clone(deep = TRUE)
func_pw_int$otu_table <- as.data.frame(round(func_pw_int$otu_table))
func_pw_int$cal_abund()

all_pathways <- rownames(func_pw$otu_table)
core_methods <- c("metagenomeSeq", "DESeq2", "ancombc2")
groups_to_test <- c("season", "site", "season_site")

# Define the standard methods (we will add GMPR/Wrench dynamically inside the loop)
base_methods <- list(
  metagenomeSeq = list(method = "metagenomeSeq"),
  linda         = list(method = "linda"),
  ALDEx2_kw     = list(method = "ALDEx2_kw"),
  ancombc2      = list(method = "ancombc2"),
  DESeq2        = list(method = "DESeq2"),
  edgeR         = list(method = "edgeR")
)

all_methods_names <- c(names(base_methods), "wilcox_GMPR", "wilcox_Wrench")

# --- Outer Loop: Iterate through each grouping variable ---

for (current_group in groups_to_test) {
  cat(sprintf("\n\n======================================================\n"))
  cat(sprintf("=== STARTING ANALYSIS FOR GROUP: %s ===\n", toupper(current_group)))
  cat(sprintf("======================================================\n"))
  
  sig_lists <- list()
  
  # 1. Run the 6 Base Methods
  for (method_name in names(base_methods)) {
    cat(sprintf("\n-> Running %s on %s...\n", method_name, current_group))
    
    current_algo <- base_methods[[method_name]]$method
    
    t_diff_res <- tryCatch({
      t_diff <- trans_diff$new(
        dataset = func_pw_int, 
        method = current_algo, 
        group = current_group, 
        p_adjust_method = "fdr"
      )
      t_diff$res_diff
    }, error = function(e) {
      cat(sprintf("   ERROR in %s: %s\n", method_name, e$message))
      return(NULL)
    })
    
    if (!is.null(t_diff_res) && nrow(t_diff_res) > 0) {
      file_path <- file.path(output_dir, sprintf("DA_%s_%s.csv", current_group, method_name))
      write.csv(t_diff_res, file_path)
      
      p_col <- if("P.adj" %in% colnames(t_diff_res)) "P.adj" else "P.value"
      sig_hits <- t_diff_res[!is.na(t_diff_res[[p_col]]) & t_diff_res[[p_col]] < 0.05, ]
      feature_col <- if("Taxa" %in% colnames(sig_hits)) "Taxa" else rownames(sig_hits)
      
      sig_lists[[method_name]] <- sig_hits[[feature_col]]
      cat(sprintf("   Done. Found %d significant pathways.\n", nrow(sig_hits)))
    } else {
      sig_lists[[method_name]] <- character(0)
    }
  }
  
  # 2. Run GMPR + Wilcoxon
  cat(sprintf("\n-> Running wilcox_GMPR on %s...\n", current_group))
  tryCatch({
    tmp_norm_gmpr <- trans_norm$new(func_pw_int)
    norm_obj_gmpr <- tmp_norm_gmpr$norm(method = "GMPR")
    # For functional data, we fake a taxonomy level so trans_diff doesn't fail
    norm_obj_gmpr$add_rownames2taxonomy("Pathway") 
    norm_obj_gmpr$cal_abund(rel = FALSE)
    
    t_diff_gmpr <- trans_diff$new(dataset = norm_obj_gmpr, method = "wilcox", group = current_group, taxa_level = "Pathway")
    
    file_path <- file.path(output_dir, sprintf("DA_%s_wilcox_GMPR.csv", current_group))
    write.csv(t_diff_gmpr$res_diff, file_path)
    
    p_col <- if("P.adj" %in% colnames(t_diff_gmpr$res_diff)) "P.adj" else "P.value"
    sig_hits <- t_diff_gmpr$res_diff[!is.na(t_diff_gmpr$res_diff[[p_col]]) & t_diff_gmpr$res_diff[[p_col]] < 0.05, ]
    feature_col <- if("Taxa" %in% colnames(sig_hits)) "Taxa" else rownames(sig_hits)
    
    sig_lists[["wilcox_GMPR"]] <- sig_hits[[feature_col]]
    cat(sprintf("   Done. Found %d significant pathways.\n", nrow(sig_hits)))
  }, error = function(e) {
    cat(sprintf("   ERROR in wilcox_GMPR: %s\n", e$message))
    sig_lists[["wilcox_GMPR"]] <- character(0)
  })

  
# Run missing Wilcoxon Methods #

output_dir <- "Output/4.PICRUSt2"
groups_to_test <- c("season", "site", "season_site")
core_methods <- c("metagenomeSeq", "DESeq2", "ancombc2")
all_methods_names <- c("metagenomeSeq", "linda", "ALDEx2_kw", "ancombc2", 
                       "DESeq2", "edgeR", "wilcox_GMPR", "wilcox_Wrench")

# --- 1. Fix the Base Object for GMPR and Wrench ---
cat("\n=== Preparing Base Object for Normalization ===\n")
func_pw_int <- func_pw$clone(deep = TRUE)
func_pw_int$otu_table <- as.data.frame(round(func_pw_int$otu_table))

# The crucial fix: Add a dummy taxonomy table so microeco's trans_norm doesn't panic
func_pw_int$tax_table <- data.frame(Pathway = rownames(func_pw_int$otu_table), 
                                    row.names = rownames(func_pw_int$otu_table))
func_pw_int$cal_abund()

all_pathways <- rownames(func_pw_int$otu_table)

# --- 2. Run ONLY the missing GMPR and Wrench methods ---
for (current_group in groups_to_test) {
  cat(sprintf("\n--- Running Missing Normalizations for %s ---\n", toupper(current_group)))
  
  # GMPR + Wilcox
  tryCatch({
    tmp_norm_gmpr <- trans_norm$new(func_pw_int)
    norm_obj_gmpr <- tmp_norm_gmpr$norm(method = "GMPR")
    norm_obj_gmpr$cal_abund(rel = FALSE)
    
    t_diff_gmpr <- trans_diff$new(dataset = norm_obj_gmpr, method = "wilcox", group = current_group, taxa_level = "Pathway")
    write.csv(t_diff_gmpr$res_diff, file.path(output_dir, sprintf("DA_%s_wilcox_GMPR.csv", current_group)))
    cat("   + wilcox_GMPR Success.\n")
  }, error = function(e) cat("   - wilcox_GMPR Failed:", e$message, "\n"))
  
  # Wrench + Wilcox
  tryCatch({
    tmp_norm_wrench <- trans_norm$new(func_pw_int)
    norm_obj_wrench <- tmp_norm_wrench$norm(method = "Wrench", condition = current_group)
    norm_obj_wrench$cal_abund(rel = FALSE)
    
    t_diff_wrench <- trans_diff$new(dataset = norm_obj_wrench, method = "wilcox", group = current_group, taxa_level = "Pathway")
    write.csv(t_diff_wrench$res_diff, file.path(output_dir, sprintf("DA_%s_wilcox_Wrench.csv", current_group)))
    cat("   + wilcox_Wrench Success.\n")
  }, error = function(e) cat("   - wilcox_Wrench Failed:", e$message, "\n"))
}


## Step 6 — FINAL CONSENSUS TABLE REBUILDER ##

output_dir <- "Output/4.PICRUSt2"
groups_to_test <- c("season", "site", "season_site")
core_methods <- c("metagenomeSeq", "DESeq2", "ancombc2")
all_methods_names <- c("metagenomeSeq", "linda", "ALDEx2_kw", "ancombc2", 
                       "DESeq2", "edgeR", "wilcox_GMPR", "wilcox_Wrench")

# Get the true list of all pathway IDs
all_pathways <- rownames(func_pw$otu_table)

cat("\n=== REBUILDING CONSENSUS TABLES (PIPED NAME FIX) ===\n")

for (current_group in groups_to_test) {
  cat(sprintf("\n-> Compiling %s...\n", current_group))
  sig_lists <- list()
  
  for (method_name in all_methods_names) {
    file_path <- file.path(output_dir, sprintf("DA_%s_%s.csv", current_group, method_name))
    
    if (file.exists(file_path)) {
      res <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Determine p-value column 
      p_cols_to_check <- c("P.adj", "padj", "P.value", "pvalue", "FDR")
      p_col <- p_cols_to_check[p_cols_to_check %in% colnames(res)][1]
      
      if (is.na(p_col)) {
        cat(sprintf("   Warning: Could not find P-value column for %s\n", method_name))
        sig_lists[[method_name]] <- character(0)
        next
      }
      
      # Filter for significance
      sig_hits <- res[!is.na(res[[p_col]]) & res[[p_col]] < 0.05, ]
      
      extracted_features <- character(0)
      
      if (nrow(sig_hits) > 0) {
        for (col in colnames(sig_hits)) {
          # THE FIX: Strip out the "|" and everything after it before checking
          cleaned_col <- gsub("\\|.*", "", as.character(sig_hits[[col]]))
          
          if (any(cleaned_col %in% all_pathways)) {
            extracted_features <- cleaned_col
            break
          }
        }
        
        # Fallback check for row names
        if (length(extracted_features) == 0) {
          cleaned_rownames <- gsub("\\|.*", "", rownames(sig_hits))
          if (any(cleaned_rownames %in% all_pathways)) {
            extracted_features <- cleaned_rownames
          }
        }
      }
      
      sig_lists[[method_name]] <- extracted_features
      cat(sprintf("   %s: Found %d significant hits.\n", method_name, length(extracted_features)))
      
    } else {
      cat(sprintf("   Warning: %s is missing for %s. Marking as 0.\n", method_name, current_group))
      sig_lists[[method_name]] <- character(0)
    }
  }
  
  # Build DataFrame
  consensus_df <- data.frame(Pathway_ID = all_pathways, stringsAsFactors = FALSE)
  if (exists("pw_desc")) consensus_df$Description <- pw_desc[consensus_df$Pathway_ID]
  
  # Add 1/0 columns
  for (method_name in all_methods_names) {
    consensus_df[[method_name]] <- ifelse(consensus_df$Pathway_ID %in% sig_lists[[method_name]], 1, 0)
  }
  
  # Summaries
  consensus_df$Total_Methods_Sig <- rowSums(consensus_df[, all_methods_names], na.rm = TRUE)
  consensus_df$Found_By_All <- ifelse(consensus_df$Total_Methods_Sig == length(all_methods_names), TRUE, FALSE)
  consensus_df$Found_By_All_But_One <- ifelse(consensus_df$Total_Methods_Sig >= (length(all_methods_names) - 1), TRUE, FALSE)
  consensus_df$Found_By_Core <- apply(consensus_df[, core_methods], 1, function(row) sum(row, na.rm=TRUE) == length(core_methods))
  
  # Sort table to put highest consensus at the top
  consensus_df <- consensus_df[order(-consensus_df$Total_Methods_Sig), ]
  
  # Save
  synth_path <- file.path(output_dir, sprintf("DA_%s_Consensus_Synthesis.csv", current_group))
  write.csv(consensus_df, synth_path, row.names = FALSE)
}



