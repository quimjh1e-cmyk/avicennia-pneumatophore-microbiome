load("Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData")
# Check how many ASVs are unresolved at genus level
tax <- amplicon_16S_microtable_rarefy$tax_table
n_total <- nrow(tax)
n_no_genus <- sum(is.na(tax$Genus) | grepl("Incertae", tax$Genus) | tax$Genus == "")
cat("ASVs without genus classification:", 9290, "/", n_total, "(", round(9290/n_total*100, 1), "%)\n")

# See unique genus names that might indicate unresolved classification
head(sort(table(tax$Genus), decreasing = TRUE), 20)

n_no_family <- sum(tax$Family == "f__")

cat("ASVs without family classification:", n_no_family, "/", n_total,"(", round(n_no_family/n_total*100, 1), "%)\n")

# Get ASVs with no genus classification
tax_no_genus <- tax[tax$Genus == "g__", ]

# Get their total abundance across all samples
otu <- amplicon_16S_microtable_rarefy$otu_table
abundance <- rowSums(otu[rownames(tax_no_genus), ])

# Top 15 most abundant unresolved ASVs with their finest classification
top_unresolved <- data.frame(
  ASV = names(sort(abundance, decreasing = TRUE))[1:15],
  Total_Reads = sort(abundance, decreasing = TRUE)[1:15],
  Family = tax_no_genus[names(sort(abundance, decreasing = TRUE))[1:15], "Family"],
  Order = tax_no_genus[names(sort(abundance, decreasing = TRUE))[1:15], "Order"]
)
print(top_unresolved)



# Get the sequences for the top 15 unresolved ASVs
top_ids <- top_unresolved$ASV
seqs <- amplicon_16S_microtable_rarefy$rep_fasta[top_ids]
# Print in FASTA format
for (i in seq_along(seqs)) {
  cat(paste0(">", names(seqs)[i]), "\n", as.character(seqs[i]), "\n")
}

library(microeco)

# Load rarefied data
#load("path/to/amplicon_16S_microtable_rarefy.RData")

# Get OTU and taxonomy tables
otu <- amplicon_16S_microtable_rarefy$otu_table
tax <- amplicon_16S_microtable_rarefy$tax_table

# Total reads across all samples
total_reads <- sum(otu)

# 1. Abundance of the 10 core genera
core_genera <- c("g__Tunicatimonas", "g__Palleronia-Pseudomaribius",
                 "g__Lewinella", "g__Gramella", "g__Ilumatobacter",
                 "g__Rubrivirga", "g__Sphingomicrobium", "g__Acuticoccus",
                 "g__Woeseia", "g__Pir4_lineage", "g__Pleurocapsa")

core_results <- data.frame(
  Genus = character(),
  N_ASVs = integer(),
  Total_Reads = integer(),
  Rel_Abundance_Pct = numeric(),
  stringsAsFactors = FALSE
)

for (g in core_genera) {
  asvs <- rownames(tax[tax$Genus == g, ])
  n_asvs <- length(asvs)
  reads <- if (n_asvs > 0) sum(otu[asvs, ]) else 0
  core_results <- rbind(core_results, data.frame(
    Genus = g, N_ASVs = n_asvs, Total_Reads = reads,
    Rel_Abundance_Pct = round(reads / total_reads * 100, 3)
  ))
}

core_results <- core_results[order(-core_results$Total_Reads), ]
cat("\n=== Core Microbiome Genera ===\n")
print(core_results, row.names = FALSE)
cat("\nTotal core abundance:", round(sum(core_results$Rel_Abundance_Pct), 2), "%\n")

# 2. Compare with the 3 BLASTed ASVs
blast_asvs <- c("ASV_3566", "ASV_10845", "ASV_1823")
blast_names <- c("Euzebyaceae (novel genus)", "Qipengyuania", "Nitriliruptor")

cat("\n=== BLASTed Unresolved ASVs ===\n")
for (i in seq_along(blast_asvs)) {
  reads <- sum(otu[blast_asvs[i], ])
  cat(sprintf("%s (%s): %d reads (%.3f%%)\n",
              blast_asvs[i], blast_names[i], reads, reads / total_reads * 100))
}

# 3. For context: total unresolved (g__) abundance
unresolved <- rownames(tax[tax$Genus == "g__", ])
unresolved_reads <- sum(otu[unresolved, ])
cat(sprintf("\n=== Unresolved (g__) total: %d reads (%.1f%%) from %d ASVs ===\n",
            unresolved_reads, unresolved_reads / total_reads * 100, length(unresolved)))

#------------------------------BLAST Pipeline analyze unassigned ASVs --------------------------------------------

library(httr)
library(xml2)
library(microeco)

# Load data
load("Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData")

otu <- amplicon_16S_microtable_rarefy$otu_table
tax <- amplicon_16S_microtable_rarefy$tax_table

# Get top 15 unresolved ASVs
unresolved_tax <- tax[tax$Genus == "g__", ]
abundance <- rowSums(otu[rownames(unresolved_tax), ])
top15 <- names(sort(abundance, decreasing = TRUE))[1:15]
seqs <- amplicon_16S_microtable_rarefy$rep_fasta[top15]

# Scoring function for hit prioritization
score_hit <- function(description, pct_identity) {
  # Must be >= 97% identity to be considered for genus assignment
  if (is.na(pct_identity) || pct_identity < 97) return(-1)
  
  desc_lower <- tolower(description)
  
  # Exclude uncultured/environmental unless nothing else
  is_uncultured <- grepl("uncultured|environmental sample|metagenome", desc_lower)
  
  # Priority scoring (higher = better)
  score <- 0
  
  if (!is_uncultured) {
    score <- score + 100  # Base bonus for cultured organism
    
    # Type strain highest priority
    if (grepl("type strain|\\(t\\)", desc_lower)) score <- score + 50
    
    # Complete genome next
    if (grepl("complete genome|complete sequence", desc_lower)) score <- score + 30
    
    # Named species (has binomial name pattern)
    if (grepl("^[A-Z][a-z]+ [a-z]+", description)) score <- score + 20
    
    # Strain identified
    if (grepl("strain ", desc_lower)) score <- score + 10
  }
  
  # Add identity as tiebreaker (0-100 scale)
  score <- score + pct_identity / 100
  
  return(score)
}

# BLAST function that retrieves top 10 hits and picks the best
blast_one_smart <- function(asv_id, sequence, retries = 3) {
  cat(sprintf("\n--- BLASTing %s ---\n", asv_id))
  
  for (attempt in 1:retries) {
    tryCatch({
      blast_url <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
      
      # Submit BLAST
      submit <- POST(blast_url, body = list(
        CMD = "Put",
        PROGRAM = "blastn",
        DATABASE = "nt",
        MEGABLAST = "on",
        QUERY = paste0(">", asv_id, "\n", as.character(sequence))
      ))
      submit_text <- content(submit, "text")
      rid <- gsub("RID = ", "", regmatches(submit_text, regexpr("RID = [A-Z0-9]+", submit_text)))
      cat(sprintf("  RID: %s\n", rid))
      
      # Poll for results
      ready <- FALSE
      for (i in 1:40) {
        Sys.sleep(15)
        check <- GET(blast_url, query = list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid))
        check_text <- content(check, "text")
        if (grepl("Status=READY", check_text)) { ready <- TRUE; break }
        if (grepl("Status=FAILED", check_text)) stop("BLAST failed")
        cat(sprintf("  Waiting... (%ds)\n", i * 15))
      }
      if (!ready) stop("BLAST timed out")
      
      # Get XML results - request 10 hits
      results_xml <- GET(blast_url, query = list(
        CMD = "Get", FORMAT_TYPE = "XML", RID = rid,
        DESCRIPTIONS = 10, ALIGNMENTS = 10
      ))
      xml_text <- content(results_xml, "text")
      doc <- read_xml(xml_text)
      
      # Parse all hits
      hits <- xml_find_all(doc, ".//Hit")
      
      if (length(hits) == 0) {
        cat("  No hits found\n")
        return(data.frame(
          ASV = asv_id, Total_Reads = sum(otu[asv_id, ]),
          SILVA_Family = tax[asv_id, "Family"],
          Best_Hit = "No hits", Identity_Pct = NA, Accession = NA,
          Selection_Reason = "No hits", All_Hits_Summary = NA,
          stringsAsFactors = FALSE
        ))
      }
      
      # Build data frame of all hits
      hit_data <- data.frame(
        description = character(), accession = character(),
        pct_identity = numeric(), score = numeric(),
        stringsAsFactors = FALSE
      )
      
      for (h in hits) {
        h_def <- xml_text(xml_find_first(h, ".//Hit_def"))
        h_acc <- xml_text(xml_find_first(h, ".//Hit_accession"))
        h_id <- as.numeric(xml_text(xml_find_first(h, ".//Hsp_identity")))
        h_len <- as.numeric(xml_text(xml_find_first(h, ".//Hsp_align-len")))
        h_pct <- round(h_id / h_len * 100, 2)
        h_score <- score_hit(h_def, h_pct)
        
        hit_data <- rbind(hit_data, data.frame(
          description = h_def, accession = h_acc,
          pct_identity = h_pct, score = h_score,
          stringsAsFactors = FALSE
        ))
      }
      
      # Sort by score (descending), then identity
      hit_data <- hit_data[order(-hit_data$score, -hit_data$pct_identity), ]
      
      # Best hit
      best <- hit_data[1, ]
      
      # Determine selection reason
      if (best$score >= 100) {
        reason <- "Cultured organism >= 97% identity"
      } else if (best$score >= 0) {
        reason <- "Uncultured only (no cultured hit >= 97%)"
      } else {
        reason <- "No hit >= 97% identity"
      }
      
      # Summarize all hits for review
      all_summary <- paste(
        sprintf("%s (%s, %.2f%%, score=%.1f)", 
                hit_data$description, hit_data$accession, 
                hit_data$pct_identity, hit_data$score),
        collapse = " | "
      )
      
      cat(sprintf("  SELECTED: %s (%s) - %.2f%% [%s]\n", 
                  best$description, best$accession, best$pct_identity, reason))
      
      # Also print runner-up if different
      if (nrow(hit_data) > 1 && hit_data$description[2] != best$description) {
        cat(sprintf("  Runner-up: %s (%.2f%%)\n", 
                    hit_data$description[2], hit_data$pct_identity[2]))
      }
      
      return(data.frame(
        ASV = asv_id,
        Total_Reads = sum(otu[asv_id, ]),
        SILVA_Family = tax[asv_id, "Family"],
        Best_Hit = best$description,
        Identity_Pct = best$pct_identity,
        Accession = best$accession,
        Selection_Reason = reason,
        All_Hits_Summary = all_summary,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat(sprintf("  Attempt %d failed: %s\n", attempt, e$message))
      if (attempt == retries) {
        return(data.frame(
          ASV = asv_id, Total_Reads = sum(otu[asv_id, ]),
          SILVA_Family = tax[asv_id, "Family"],
          Best_Hit = "FAILED", Identity_Pct = NA, Accession = NA,
          Selection_Reason = "BLAST failed", All_Hits_Summary = NA,
          stringsAsFactors = FALSE
        ))
      }
      Sys.sleep(30)
    })
  }
}

# Run pipeline
cat("Starting smart BLAST pipeline for 15 ASVs...\n")
cat("Estimated time: 20-30 minutes\n\n")

results_list <- list()
for (i in seq_along(top15)) {
  cat(sprintf("\n===== [%d/15] =====", i))
  results_list[[i]] <- blast_one_smart(top15[i], seqs[top15[i]])
  Sys.sleep(5)
}

blast_results <- do.call(rbind, results_list)

# Save full results
write.csv(blast_results,
          "Output/1.Amplicon/taxonomy/blast_top15_smart_results.csv",
          row.names = FALSE)

# Print clean summary
cat("\n\n========== FINAL RESULTS ==========\n\n")
summary_df <- blast_results[, c("ASV", "Total_Reads", "SILVA_Family", "Best_Hit", "Identity_Pct", "Accession", "Selection_Reason")]
print(summary_df, row.names = FALSE, right = FALSE)

cat("\nResults saved to blast_top15_smart_results.csv\n")

# ============================================================
# FIGURE 3 — Rivularia PCC-7116 relative abundance
# Computed from unrarefied microtable: all 16 ASVs classified
# as g__Rivularia_PCC-7116 in SILVA taxonomy.
# Rivularia is completely absent from Cola de Ballena (both
# seasons) but present at Estero Zacatecas, supporting the
# mismatch between PICRUSt2 functional predictions and 16S
# taxonomy discussed in the thesis text.
# Output: figures/rivularia_abundance.pdf / .png
# ============================================================

# ── F3-1. Compute relative abundance from microtable ────────
load("Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable.RData")

riv_asvs <- rownames(amplicon_16S_microtable$tax_table)[
  grep("Rivularia", amplicon_16S_microtable$tax_table$Genus)]
riv_counts <- colSums(amplicon_16S_microtable$otu_table[riv_asvs, , drop = FALSE])
sample_totals <- colSums(amplicon_16S_microtable$otu_table)
riv_rel <- (riv_counts / sample_totals) * 100

sample_order <- c(paste0("FCBAG", 1:6), paste0("FEZAG", 1:6),
                  paste0("SCBAG", 1:6), paste0("SEZAG", 1:6))

riv_df <- data.frame(
  sample = sample_order,
  rel_ab = riv_rel[sample_order],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Season = ifelse(grepl("^F", sample), "February", "September"),
    Site   = ifelse(grepl("CBA", sample), "Cola de Ballena", "Estero Zacatecas"),
    group  = paste(Season, Site, sep = "\n")
  )

message(sprintf("Rivularia ASVs found: %d", length(riv_asvs)))

group_levels <- c(
  "February\nCola de Ballena",
  "February\nEstero Zacatecas",
  "September\nCola de Ballena",
  "September\nEstero Zacatecas"
)
riv_df$group <- factor(riv_df$group, levels = group_levels)


# ── F3-2. Build plot (bar chart, all 24 samples) ────────────
# Sample order matching the heatmap: FCBAG1-6, FEZAG1-6, SCBAG1-6, SEZAG1-6
riv_df$sample <- factor(riv_df$sample,
  levels = c(paste0("FCBAG", 1:6), paste0("FEZAG", 1:6),
             paste0("SCBAG", 1:6), paste0("SEZAG", 1:6)))

riv_df$Group <- paste(riv_df$Season, riv_df$Site, sep = " — ")
group_pal <- c(
  "February — Cola de Ballena"      = "#56B4E9",
  "February — Estero Zacatecas"     = "#CC79A7",
  "September — Cola de Ballena"     = "#009E73",
  "September — Estero Zacatecas"    = "#E69F00"
)

riv_plot <- ggplot(riv_df, aes(x = sample, y = rel_ab, fill = Group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = group_pal) +
  labs(
    x    = NULL,
    y    = expression(italic("Rivularia")~"PCC-7116 relative abundance (%)"),
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x        = element_text(size = 7.5, angle = 45, hjust = 1, colour = "black"),
    axis.text.y        = element_text(size = 8, colour = "black"),
    axis.title.y       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.35, "cm"),
    legend.position    = "top",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ── F3-3. Save ───────────────────────────────────────────────
ggsave("Output/1.Amplicon/taxonomy/rivularia_abundance.pdf", riv_plot,
       width = 8, height = 4, device = cairo_pdf)
message("Saved: figures/rivularia_abundance.pdf")

ggsave("Output/1.Amplicon/taxonomy/rivularia_abundance.png", riv_plot,
       width = 8, height = 4, dpi = 300)
message("Saved: figures/rivularia_abundance.png")

# ── F3b. Stacked bar chart: per-ASV breakdown ───────────────
# Shows which of the 16 Rivularia ASVs dominate each sample.
# Top ASVs by total abundance get individual colours; the rest
# are pooled as "Other".

riv_mat <- amplicon_16S_microtable$otu_table[riv_asvs, sample_order, drop = FALSE]
riv_rel_mat <- sweep(riv_mat, 2, sample_totals[sample_order], "/") * 100

asv_totals <- rowSums(riv_mat)
top_n <- 6
top_asvs <- names(sort(asv_totals, decreasing = TRUE))[1:top_n]

riv_long <- as.data.frame(as.table(as.matrix(riv_rel_mat)))
colnames(riv_long) <- c("ASV", "sample", "rel_ab")
riv_long$ASV <- as.character(riv_long$ASV)
riv_long$sample <- as.character(riv_long$sample)
riv_long$ASV_label <- ifelse(riv_long$ASV %in% top_asvs, riv_long$ASV, "Other")

riv_stacked <- riv_long %>%
  group_by(sample, ASV_label) %>%
  summarise(rel_ab = sum(rel_ab), .groups = "drop")

asv_levels <- c(sort(top_asvs), "Other")
riv_stacked$ASV_label <- factor(riv_stacked$ASV_label, levels = asv_levels)
riv_stacked$sample <- factor(riv_stacked$sample, levels = sample_order)

asv_pal <- c(
  setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628"),
    sort(top_asvs)),
  Other = "grey70"
)

riv_stacked_plot <- ggplot(riv_stacked, aes(x = sample, y = rel_ab, fill = ASV_label)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = asv_pal) +
  labs(
    x    = NULL,
    y    = expression(italic("Rivularia")~"PCC-7116 relative abundance (%)"),
    fill = "ASV"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x        = element_text(size = 7.5, angle = 45, hjust = 1, colour = "black"),
    axis.text.y        = element_text(size = 8, colour = "black"),
    axis.title.y       = element_text(size = 9),
    legend.text        = element_text(size = 7.5),
    legend.key.size    = unit(0.3, "cm"),
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("Output/1.Amplicon/taxonomy/rivularia_abundance_by_ASV.pdf", riv_stacked_plot,
       width = 9, height = 4.5, device = cairo_pdf)
message("Saved: figures/rivularia_abundance_by_ASV.pdf")

ggsave("Output/1.Amplicon/taxonomy/rivularia_abundance_by_ASV.png", riv_stacked_plot,
       width = 9, height = 4.5, dpi = 300)
message("Saved: figures/rivularia_abundance_by_ASV.png")