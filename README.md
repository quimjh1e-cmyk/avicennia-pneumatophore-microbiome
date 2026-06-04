# Taxonomic Diversity and Metabolic Potential of Microbial Communities on *Avicennia germinans* Pneumatophores

Master's Thesis (TFM) from the Universitat Oberta de Catalunya (UOC), 2026

## Overview

This repository contains the bioinformatic analysis code and processed data for a 16S rRNA amplicon sequencing study of prokaryotic communities colonizing pneumatophores (aerial roots) of the black mangrove *Avicennia germinans* in two sites of Baja California Sur, Mexico.

**Study design:**
- 2 sites: Cola de Ballena (CB) and Estero Zacatecas (EZ)
- 2 seasons: February 2024 (cold, dry) and September 2024 (warm, rainy)
- 6 trees per site per season = 24 samples
- 16S rRNA V3-V4 region, paired-end Illumina sequencing

## Repository Structure

```
├── data/
│   ├── metadata/          Sample metadata and QIIME2 manifest files
│   └── abiotic/           Environmental/physicochemical measurements
├── analysis/
│   ├── microeco/          R analysis pipeline (microeco package)
│   │   ├── Rcode/        Step-by-step R scripts (Steps 1-22)
│   │   ├── Input/        Processed ASV table, taxonomy, tree
│   │   └── Output/       Result figures, statistical tables, and taxonomy/BLAST verification
│   └── picrust2/          Functional prediction analysis
│       ├── run_picrust2_pipeline.sh   Reusable pipeline script
│       ├── scripts/       R scripts for heatmaps and pathway-taxa analysis
│       ├── input/         ASV table and sequences for PICRUSt2
│       ├── output/        PICRUSt2 predicted pathway/KO/EC abundances
│       └── results/       Figures and pathway-taxa contribution tables
└── scripts/               Utility shell scripts for DA searches
```

## Analysis Pipeline

### 1. Sequence Processing (QIIME2 v2026.1.0)

Raw FASTQ files were processed in QIIME2:
- Primer trimming with cutadapt
- Denoising and ASV generation with DADA2
- Taxonomy assignment with a pre-trained SILVA 138.2 V3-V4 classifier
- Chloroplast and mitochondria filtering
- Phylogenetic tree construction (MAFFT + FastTree)

### 2. Community Analysis (R — microeco)

The `analysis/microeco/Rcode/` directory contains numbered scripts that should be run sequentially:

| Script | Analysis |
|--------|----------|
| Step 1 | Install packages |
| Step 2 | Import QIIME2 data into microtable |
| Step 3 | Preprocess microtable |
| Step 4 | Normalization (multiple methods) |
| Steps 5-6 | Core microbiome identification |
| Steps 7-9 | Alpha diversity and environmental correlations |
| Steps 10-11 | Beta diversity (PCoA, PERMANOVA, ANOSIM) |
| Steps 12-18 | Differential abundance (DESeq2, ANCOM-BC2, LEfSe, metagenomeSeq, etc.) |
| Steps 19-22 | Machine learning classification and regression |

Set your R working directory to `analysis/microeco/` before running.

### 3. Functional Prediction (PICRUSt2)

Metabolic potential predicted from 16S ASV data using PICRUSt2, with focus on nitrogen cycling (N-fixation, denitrification), carbon cycling (CODH pathway), and phosphorus and sulfur metabolism.

The pipeline can be reproduced with the included script:

```bash
conda activate picrust2
./analysis/picrust2/run_picrust2_pipeline.sh \
    -i analysis/picrust2/input \
    -o analysis/picrust2/output \
    -t 4 -n 2
```

This runs the full PICRUSt2 workflow: BIOM conversion, phylogenetic placement, gene family and pathway prediction, NSTI quality control, and description annotation. R scripts in `analysis/picrust2/scripts/` handle downstream visualization (heatmaps, pathway-taxa contributions).

### 4. Differential Abundance Search Scripts

```bash
# Search for a taxon across all DA methods (season comparison)
./scripts/search_da_season.sh Rivularia

# Search for a taxon across all DA methods (site comparison)
./scripts/search_da_site.sh Prolixibacteraceae
```

## Raw Sequence Data

Raw FASTQ files (48 paired-end files, ~957 MB) are deposited in NCBI SRA:
- BioProject: *[accession pending]*
- SRA: *[accession pending]*

## Software Requirements

- **QIIME2** v2026.1.0 (sequence processing)
- **R** >= 4.3 with packages:
  - microeco, file2meco
  - vegan, phyloseq
  - DESeq2, ANCOMBC, metagenomeSeq, ALDEx2, edgeR, LinDA, MaAsLin2, lefser
  - ggplot2, ComplexHeatmap
- **PICRUSt2** v2.5+ (functional prediction)

## Citation

Heintze, Q. (2025). Spatiotemporal variation and anthropogenic influence on epiphytic prokaryotic communities of *Avicennia germinans* pneumatophores. Master's Thesis, Universitat Oberta de Catalunya.

## License

Code in this repository is available under the MIT License. The microeco protocol scripts in `analysis/microeco/` are adapted from [Liu et al. (2025)](https://doi.org/10.1038/s41596-024-01133-1) — see `analysis/microeco/LICENSE` for their original license.
