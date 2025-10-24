# MT-2 and MT-3 ISS Microbiome Analysis

R scripts for analyzing microbiome data from the International Space Station (ISS) Missions MT-2 and MT-3, including differential abundance analysis, antibiotic resistance genes (ARGs), virulence factors (VFs), lateral gene transfer (LGT), and phenotypic characterization.

## Contents

### Analysis Scripts

1. **01_differential_abundance_analysis.R**
   - ALDEx2 differential abundance testing
   - Mission comparisons (MT-2 vs MT-3)
   - PMA vs untreated comparisons
   - Surface-specific analyses
   - Flight and collection day comparisons

2. **02_ARG_VF_LGT_analysis.R**
   - Antibiotic resistance gene (ARG) abundance analysis
   - Virulence factor (VF) quantification
   - Lateral gene transfer (LGT) event detection
   - ARG mechanism of action analysis
   - Statistical comparisons across missions, flights, surfaces, and collection days

3. **03_visualization_plots.R**
   - Combined ARG/VF/LGT visualizations
   - ARG mechanism bubble plots
   - Antibiotic resistance heatmaps
   - Virulence survival curves
   - Predicted phenotype bar charts

4. **04_taxa_abundance_skin_qPCR.R**
   - Species abundance summaries
   - Skin-associated microbiome analysis
   - qPCR validation analysis
   - Species of concern tracking
   - Phyloseq object creation

## Dependencies

### Required R Packages

```r
# Differential abundance
install.packages("ALDEx2")

# Microbiome analysis
install.packages("phyloseq")
install.packages("mctoolsr")
install.packages("vegan")

# Data manipulation
install.packages("dplyr")
install.packages("tidyr")

# Visualization
install.packages("ggplot2")
install.packages("patchwork")
install.packages("pheatmap")
install.packages("viridis")
install.packages("RColorBrewer")
install.packages("forcats")
install.packages("hrbrthemes")

# Statistics
install.packages("rstatix")
install.packages("ggpubr")
```

## Required Input Files

### Microbiome Data
- `physeq_feature_table_mt_unrar.txt` - Unrarefied OTU/ASV table
- `physeq_taxonomy_mt_unrar.txt` - Taxonomy assignments
- `MT-2_MT-3_combined_metadata_31Jan.txt` - Sample metadata
- `all_species_abs_abundance_table_MSH.txt` - Species-level abundance table

### ARG/VF/LGT Data
- `all_mech_MT2_MSH_13Nov2024.txt` - MT-2 ARG mechanism data
- `all_mech_MT3_MSH_13Nov2024.txt` - MT-3 ARG mechanism data

### Phenotype Data
- `Isolate_AMR_Heatmap.txt` - Isolate antibiotic resistance profiles
- `most_virulent_all_missions.txt` - Virulence assay time-series data
- `percent_pheno.txt` - Predicted phenotype proportions

### qPCR Data
- `MT3_qPCR_all_unmassaged.txt` - qPCR quantification data

### Rarefied Data (optional, will be generated)
- `MT2_MT3_SBG_rar.txt` - Rarefied abundance table
- `physeq_feature_table_MT2_MT3_bacteria_rar.txt` - Rarefied phyloseq table
- `physeq_taxonomy_MT2_MT3_bacteria_rar.txt` - Rarefied taxonomy

## Usage

### 1. Set Up Working Directory

```r
# Set working directory
setwd("/path/to/your/data")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
```

### 2. Run Analysis Scripts

Run scripts in order:

```r
# Differential abundance analysis
source("01_differential_abundance_analysis.R")

# ARG/VF/LGT analysis
source("02_ARG_VF_LGT_analysis.R")

# Generate visualizations
source("03_visualization_plots.R")

# Taxa abundance and qPCR
source("04_taxa_abundance_skin_qPCR.R")
```

### 3. Expected Outputs

#### Results Directory (`results/`)
- `wilcox_ARG_mission_noPMA.csv` - ARG mission comparisons
- `wilcox_VF_surface_noPMA.csv` - VF surface comparisons
- `wilcox_LGT_surface_noPMA.csv` - LGT surface comparisons
- `wilcox_ARG_mechanism_significant.csv` - Significant ARG mechanisms
- `qpcr_*.txt` - qPCR statistical results
- `significant_taxa_*.csv` - Differentially abundant taxa
- `aldex2_results_*.tsv` - Full ALDEx2 results

#### Figures Directory (`figures/`)
- `combined_ARG_VF_LGT.pdf` - Combined ARG/VF/LGT plots
- `ARG_mechanism_bubbleplot.pdf` - ARG mechanism visualization
- `ARG_resistance_heatmap.pdf` - Antibiotic resistance heatmap
- `virulence_survival_curves.pdf` - Virulence assay results
- `predicted_phenotype.pdf` - Phenotype prediction accuracy
- `skin_associated_taxa.pdf` - Skin microbiome analysis

#### Data Files
- `mt_physeq_rarefied.rds` - Phyloseq object (R binary)
- `sum_spp*.txt` - Species abundance summaries
- `MT2_MT3_SBG_rar.txt` - Rarefied count table

## Key Analyses

### Differential Abundance (Script 01)
- **Mission comparison**: Identifies taxa enriched in MT-2 vs MT-3
- **PMA treatment**: Evaluates impact of viability treatment
- **Surface analysis**: Pairwise comparisons against Cupola baseline
- **Temporal trends**: Flight and collection day effects

### ARG/VF/LGT (Script 02)
- **ARG abundance**: Normalized counts per million (CPM)
- **VF abundance**: Virulence factor quantification
- **LGT events**: Lateral gene transfer detection
- **Mechanism analysis**: ARG resistance mechanisms by surface/mission

### Visualization (Script 03)
- **Publication-quality figures**: PDF format, ready for manuscripts
- **Integrated plots**: Combined ARG/VF/LGT panels
- **Heatmaps**: Hierarchical clustering of resistance profiles
- **Time-series**: Virulence survival curves

### Taxa & qPCR (Script 04)
- **Rarefaction**: Normalizes sequencing depth to 10,000 reads
- **Skin microbiome**: 27 skin-associated taxa tracked
- **qPCR validation**: 16S and ITS quantification
- **Phyloseq integration**: Compatible with other microbiome tools

## Statistical Methods

- **ALDEx2**: Compositional data analysis with CLR transformation
- **Wilcoxon rank-sum**: Non-parametric pairwise comparisons
- **Kruskal-Wallis**: Multi-group comparisons
- **Benjamini-Hochberg**: FDR correction for multiple testing
- **Effect sizes**: Cohen's d for practical significance

## Notes

- **Rarefaction depth**: 10,000 reads per sample
- **Significance threshold**: p < 0.05 (FDR-adjusted)
- **Effect size threshold**: |effect| > 0.5 for biological relevance
- **Missing data**: Samples with NA surfaces are excluded
- **Flight 8**: Excluded from all analyses
- **PMA treatment**: "Yes" = viable cells, "No" = total DNA

## Citation

If using these scripts, please cite:

[Your publication details here]

## Contact

[Your contact information]

## License

[License information]
