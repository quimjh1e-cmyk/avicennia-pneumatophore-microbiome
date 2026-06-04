##
## Perform core taxa analysis for each group (season x site)
##
## Core taxa are identified at the ASV level (freq >= 50%, mean rel. abund >= 0.05%).
## Taxonomic-level visualization (Family, Genus) is handled in Step 6.
##


######################################################
# load packages
library(microeco)
library(magrittr)
library(ggplot2)
######################################################
# create an output directory if it does not exist
output_dir <- "./Output/1.Amplicon/Stage3_coremicrobiome"
if(! dir.exists(output_dir)){
	dir.create(output_dir, recursive = TRUE)
}
# load rarefied data
input_path <- "./Output/1.Amplicon/Stage2_amplicon_microtable/amplicon_16S_microtable_rarefy.RData"
# first check whether data path exists
if(! file.exists(input_path)){
	stop("Please first run all the scripts in Stage2 !")
}
load(input_path)
######################################################


# use temporary name tmp_microtable for convenience
tmp_microtable <- clone(amplicon_16S_microtable_rarefy)


# define core taxa: occurrence frequency 50%; mean relative abundance 0.05%
freq <- 0.5
abund <- 0.0005

# select core taxa for each group (season x site combination)
# February - Cola de Ballena
FCBA <- clone(tmp_microtable)
FCBA$sample_table %<>% .[.$season == "February" & .$site == "ColaDeBallena", ]
FCBA$tidy_dataset()
FCBA$filter_taxa(rel_abund = abund, freq = freq)

# February - Estero Zacatecas
FEZA <- clone(tmp_microtable)
FEZA$sample_table %<>% .[.$season == "February" & .$site == "EsteroZacatecas", ]
FEZA$tidy_dataset()
FEZA$filter_taxa(rel_abund = abund, freq = freq)

# September - Cola de Ballena
SCBA <- clone(tmp_microtable)
SCBA$sample_table %<>% .[.$season == "September" & .$site == "ColaDeBallena", ]
SCBA$tidy_dataset()
SCBA$filter_taxa(rel_abund = abund, freq = freq)

# September - Estero Zacatecas
SEZA <- clone(tmp_microtable)
SEZA$sample_table %<>% .[.$season == "September" & .$site == "EsteroZacatecas", ]
SEZA$tidy_dataset()
SEZA$filter_taxa(rel_abund = abund, freq = freq)

# merge the core ASV names into a data.frame object and save it into a file
res <- rbind(
	data.frame(group = "February_ColaDeBallena",  ASV = FCBA$taxa_names()),
	data.frame(group = "February_EsteroZacatecas", ASV = FEZA$taxa_names()),
	data.frame(group = "September_ColaDeBallena",  ASV = SCBA$taxa_names()),
	data.frame(group = "September_EsteroZacatecas", ASV = SEZA$taxa_names())
)
# write the data.frame object containing core ASV names in each group to the output directory
write.csv(res, file.path(output_dir, "Coretaxa_calc_groups.csv"))
# Description of the result file:
#	The "group" column represents the season x site combination, and the "ASV" column contains
#	the names of core ASVs identified in each group (freq >= 50%, mean rel. abund >= 0.05%).


# also save each microtable object to local files with .RData format
save(FCBA, file = file.path(output_dir, "Coretaxa_February_ColaDeBallena.RData"),  compress = TRUE)
save(FEZA, file = file.path(output_dir, "Coretaxa_February_EsteroZacatecas.RData"), compress = TRUE)
save(SCBA, file = file.path(output_dir, "Coretaxa_September_ColaDeBallena.RData"),  compress = TRUE)
save(SEZA, file = file.path(output_dir, "Coretaxa_September_EsteroZacatecas.RData"), compress = TRUE)
