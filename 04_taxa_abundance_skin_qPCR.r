################################################################################
# Taxa Abundance, Skin-Associated Microbiome, and qPCR Analysis
# MT-2 and MT-3 Microbiome Data
################################################################################

# Load required libraries
library(mctoolsr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(viridis)
library(rstatix)

################################################################################
# 1. DATA IMPORT AND RAREFACTION
################################################################################

# Import unrarefied dataset
mt <- load_taxa_table("all_species_abs_abundance_table_MSH.txt", 
                      "MT-2_MT-3_combined_metadata_31Jan.txt")

cat("Unrarefied sample sizes:\n")
print(sort(colSums(mt$data_loaded)))

# Rarefy to 10,000 reads
mt_rar <- single_rarefy(mt, 10000)

cat("\nRarefied sample sizes:\n")
print(sort(colSums(mt_rar$data_loaded)))

# Export rarefied dataset
export_taxa_table(mt_rar, "MT2_MT3_SBG_rar.txt")

# Filter out Flight 8
mt_rar <- filter_data(mt_rar, 'Flight', filter_vals = '8')

# Separate by PMA treatment
mt_rar_noPMA <- filter_data(mt_rar, 'PMA_Treated', keep_vals = 'No')
mt_rar_PMA <- filter_data(mt_rar, 'PMA_Treated', keep_vals = 'Yes')

################################################################################
# 2. SPECIES ABUNDANCE SUMMARIES
################################################################################

cat("\n=== Calculating Species Abundance ===\n")

# All samples
tax_sum_spp <- summarize_taxonomy(mt_rar, level = 7, report_higher_tax = FALSE)
sum_spp <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_spp, 
                                       'PMA_Treated', 
                                       metadata_map = mt_rar$map_loaded, 
                                       smry_fun = mean)
write.table(sum_spp, file = "sum_spp.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# Untreated samples
tax_sum_spp_noPMA <- summarize_taxonomy(mt_rar_noPMA, level = 7, 
                                         report_higher_tax = FALSE)
sum_spp_noPMA <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_spp_noPMA, 
                                             'PMA_Treated', 
                                             metadata_map = mt_rar_noPMA$map_loaded, 
                                             smry_fun = mean)
write.table(sum_spp_noPMA, file = "sum_spp_noPMA.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# PMA-treated samples
tax_sum_spp_PMA <- summarize_taxonomy(mt_rar_PMA, level = 7, 
                                       report_higher_tax = FALSE)
sum_spp_PMA <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_spp_PMA, 
                                           'PMA_Treated', 
                                           metadata_map = mt_rar_PMA$map_loaded, 
                                           smry_fun = mean)
write.table(sum_spp_PMA, file = "sum_spp_PMA.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

# PMA-treated by mission
mt_rar_PMA_MT2 <- filter_data(mt_rar_PMA, 'Mission', keep_vals = 'MT-2')
tax_sum_spp_PMA_MT2 <- summarize_taxonomy(mt_rar_PMA_MT2, level = 7, 
                                           report_higher_tax = FALSE)
sum_spp_PMA_MT2 <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_spp_PMA_MT2, 
                                               'PMA_Treated', 
                                               metadata_map = mt_rar_PMA_MT2$map_loaded, 
                                               smry_fun = mean)
write.table(sum_spp_PMA_MT2, file = "sum_spp_PMA_MT2.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

mt_rar_PMA_MT3 <- filter_data(mt_rar_PMA, 'Mission', keep_vals = 'MT-3')
tax_sum_spp_PMA_MT3 <- summarize_taxonomy(mt_rar_PMA_MT3, level = 7, 
                                           report_higher_tax = FALSE)
sum_spp_PMA_MT3 <- taxa_summary_by_sample_type(taxa_smry_df = tax_sum_spp_PMA_MT3, 
                                               'PMA_Treated', 
                                               metadata_map = mt_rar_PMA_MT3$map_loaded, 
                                               smry_fun = mean)
write.table(sum_spp_PMA_MT3, file = "sum_spp_PMA_MT3.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

################################################################################
# 3. SKIN-ASSOCIATED MICROBIOME ANALYSIS
################################################################################

cat("\n=== Analyzing Skin-Associated Taxa ===\n")

# Define skin-associated taxa
skin_taxa <- c(
  's__Acinetobacter_johnsonii',
  's__Acinetobacter_lwoffii',
  's__Corynebacterium_accolens',
  's__Corynebacterium_amycolatum',
  's__Corynebacterium_jeikeium',
  's__Corynebacterium_kroppenstedtii',
  's__Corynebacterium_mucifaciens',
  's__Corynebacterium_striatum',
  's__Corynebacterium_tuberculostearicum',
  's__Cutibacterium_acnes',
  's__Cutibacterium_avidum',
  's__Cutibacterium_granulosum',
  's__Cutibacterium_modestum',
  's__Cutibacterium_namnetense',
  's__Finegoldia_magna',
  's__Lawsonella_clevelandensis',
  's__Micrococcus_luteus',
  'g__Peptoniphilus',
  's__Staphylococcus_aureus',
  's__Staphylococcus_capitis',
  's__Staphylococcus_epidermidis',
  's__Staphylococcus_haemolyticus',
  's__Staphylococcus_hominis',
  's__Staphylococcus_lugdunensis',
  's__Staphylococcus_simulans',
  's__Staphylococcus_warneri',
  's__Staphylococcus_xylosus'
)

## UNTREATED SAMPLES ##

skin_nopma <- filter_taxa_from_input(mt_rar_noPMA, taxa_to_keep = skin_taxa)

skin.tot_nopma <- as.data.frame(colSums(skin_nopma$data_loaded))
names(skin.tot_nopma)[1] <- "Total_Abundance"

# Convert to percent
skin.per_nopma <- (skin.tot_nopma / 10000) * 100

# Merge with metadata
skin.per.map_nopma <- merge(skin.per_nopma, mt_rar_noPMA$map_loaded, 
                            by = 0, stringsAsFactors = FALSE)

# Remove NA surfaces
skin.per.map_clean_noPMA <- skin.per.map_nopma %>%
  filter(!is.na(Surface), Surface != "NA", Mission != "MT-2")

# Statistical testing
cat("\nSkin-associated taxa - Mission comparison (untreated):\n")
print(aggregate(Total_Abundance ~ Mission, skin.per.map_nopma, mean))
print(wilcox.test(Total_Abundance ~ Mission, data = skin.per.map_nopma))

cat("\nSkin-associated taxa - Surface comparison (MT-3, untreated):\n")
wilcox_skin_surface_noPMA <- skin.per.map_clean_noPMA %>%
  pairwise_wilcox_test(Total_Abundance ~ Surface, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001  ~ "***",
    p.adj <= 0.01   ~ "**",
    p.adj <= 0.05   ~ "*",
    TRUE            ~ "ns"
  ))
print(wilcox_skin_surface_noPMA)

# Plot
p1 <- ggplot(skin.per.map_clean_noPMA, 
             aes(x = Surface, y = Total_Abundance, fill = PMA_Treated)) +
  geom_boxplot() +
  labs(x = "Surface", y = "Percent relative abundance", fill = "PMA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("royalblue3", "firebrick3"))

## PMA-TREATED SAMPLES ##

mt_rar_PMA_filtered <- filter_data(mt_rar_PMA, 'Mission', keep_vals = 'MT-3')
mt_rar_PMA_filtered <- filter_data(mt_rar_PMA_filtered, 'Surface_ID', 
                                   filter_vals = c('Environmental Control', 
                                                  'Kit Control'))

skin_PMA <- filter_taxa_from_input(mt_rar_PMA_filtered, taxa_to_keep = skin_taxa)

skin.tot_PMA <- as.data.frame(colSums(skin_PMA$data_loaded))
names(skin.tot_PMA)[1] <- "Total_Abundance"

skin.per_PMA <- (skin.tot_PMA / 10000) * 100

skin.per.map_PMA <- merge(skin.per_PMA, mt_rar_PMA_filtered$map_loaded, 
                          by = 0, stringsAsFactors = FALSE)

cat("\nSkin-associated taxa - Surface comparison (MT-3, PMA-treated):\n")
wilcox_results_skin_PMA <- skin.per.map_PMA %>%
  pairwise_wilcox_test(Total_Abundance ~ Surface, p.adjust.method = "BH")
print(wilcox_results_skin_PMA)

# Plot
p2 <- ggplot(skin.per.map_PMA, aes(x = Surface, y = Total_Abundance)) +
  geom_boxplot(fill = "burlywood2") +
  labs(x = "Surface", y = "Percent relative abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combined plot
skin_combined <- p1 / p2
ggsave("figures/skin_associated_taxa.pdf", skin_combined, 
       width = 10, height = 8)

################################################################################
# 4. qPCR ANALYSIS
################################################################################

cat("\n=== qPCR Analysis ===\n")

# Load qPCR data
qpcr <- read.delim(file = "MT3_qPCR_all_unmassaged.txt", header = TRUE)

# Separate by treatment
data_nopma <- qpcr %>% filter(Treatment == "NO_PMA")
data_pma <- qpcr %>% filter(Treatment == "PMA")

# Test differences between assays
cat("\nKruskal-Wallis test - Assay comparison (untreated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Assay, data = data_nopma))

cat("\nKruskal-Wallis test - Assay comparison (PMA-treated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Assay, data = data_pma))

# Separate by assay type
data_16S <- qpcr %>% filter(Assay == "16S")
data_ITS <- qpcr %>% filter(Assay == "ITS")

data_16S_nopma <- data_16S %>% filter(Treatment == "NO_PMA")
data_16S_pma <- data_16S %>% filter(Treatment == "PMA")
data_ITS_nopma <- data_ITS %>% filter(Treatment == "NO_PMA")
data_ITS_pma <- data_ITS %>% filter(Treatment == "PMA")

# Location comparisons
cat("\n16S - Location comparison (untreated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Location, data = data_16S_nopma))

cat("\n16S - Location comparison (PMA-treated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Location, data = data_16S_pma))

cat("\nITS - Location comparison (untreated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Location, data = data_ITS_nopma))

cat("\nITS - Location comparison (PMA-treated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Location, data = data_ITS_pma))

# Pairwise comparisons
data_16S_wilcoxon_nopma <- pairwise.wilcox.test(
  data_16S_nopma$Log10_ng_per_wipe,
  data_16S_nopma$Location,
  p.adjust.method = "BH"
)

data_16S_wilcoxon_pma <- pairwise.wilcox.test(
  data_16S_pma$Log10_ng_per_wipe,
  data_16S_pma$Location,
  p.adjust.method = "BH"
)

data_ITS_wilcoxon_nopma <- pairwise.wilcox.test(
  data_ITS_nopma$Log10_ng_per_wipe,
  data_ITS_nopma$Location,
  p.adjust.method = "BH"
)

data_ITS_wilcoxon_pma <- pairwise.wilcox.test(
  data_ITS_pma$Log10_ng_per_wipe,
  data_ITS_pma$Location,
  p.adjust.method = "BH"
)

# Export results
write.table(data_16S_wilcoxon_nopma$p.value, 
            file = "results/qpcr_16S_wilcoxon_nopma.txt", 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

write.table(data_16S_wilcoxon_pma$p.value, 
            file = "results/qpcr_16S_wilcoxon_pma.txt", 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

write.table(data_ITS_wilcoxon_nopma$p.value, 
            file = "results/qpcr_ITS_wilcoxon_nopma.txt", 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

write.table(data_ITS_wilcoxon_pma$p.value, 
            file = "results/qpcr_ITS_wilcoxon_pma.txt", 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Flight comparisons
cat("\n16S - Flight comparison (untreated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Flight, data = data_16S_nopma))

cat("\n16S - Flight comparison (PMA-treated):\n")
print(kruskal.test(Log10_ng_per_wipe ~ Flight, data = data_16S_pma))

data_16S_wilcoxon_nopma_flight <- pairwise.wilcox.test(
  data_16S_nopma$Log10_ng_per_wipe,
  data_16S_nopma$Flight,
  p.adjust.method = "BH"
)

data_16S_wilcoxon_pma_flight <- pairwise.wilcox.test(
  data_16S_pma$Log10_ng_per_wipe,
  data_16S_pma$Flight,
  p.adjust.method = "BH"
)

################################################################################
# 5. SPECIES OF CONCERN ABUNDANCE
################################################################################

cat("\n=== Species of Concern Analysis ===\n")

# Filter for untreated samples without controls
mt_rar_noPMA_filtered <- filter_data(mt_rar_noPMA, 'Surface_ID', 
                                     filter_vals = c('Environmental Control', 
                                                    'Kit Control'))

# Example: Enterococcus faecalis
species_of_concern <- filter_taxa_from_input(
  mt_rar_noPMA_filtered, 
  taxa_to_keep = 's__Enterococcus_faecalis'
)

species.tot <- as.data.frame(colSums(species_of_concern$data_loaded))
names(species.tot)[1] <- "Total_Abundance"

species.per <- (species.tot / 10000) * 100

species.per.map <- merge(species.per, mt_rar_noPMA_filtered$map_loaded, 
                         by = 0, stringsAsFactors = FALSE)

cat("\nEnterococcus faecalis abundance by Mission:\n")
print(aggregate(Total_Abundance ~ Mission, species.per.map, mean))

################################################################################
# 6. PHYLOSEQ OBJECT CREATION (for additional analyses)
################################################################################

# Read metadata
mt_metadata <- read.delim("MT-2_MT-3_combined_metadata_31Jan.txt", 
                         header = TRUE, sep = '\t', row.names = 1)

# Read taxonomy
mt_taxonomy_rar <- read.delim("physeq_taxonomy_MT2_MT3_bacteria_rar.txt", 
                              header = TRUE, sep = '\t')
row.names(mt_taxonomy_rar) <- mt_taxonomy_rar$OTU_ID

# Read OTU table
mt_sbg_rar <- read.delim("physeq_feature_table_MT2_MT3_bacteria_rar.txt", 
                         header = TRUE, sep = '\t', check.names = FALSE)
if("OTU_ID" %in% colnames(mt_sbg_rar)) {
  rownames(mt_sbg_rar) <- mt_sbg_rar$OTU_ID
  mt_sbg_rar <- mt_sbg_rar[, !colnames(mt_sbg_rar) %in% c("OTU_ID")]
}

# Convert to phyloseq objects
mt_sbg_table_rar <- otu_table(as.matrix(mt_sbg_rar), taxa_are_rows = TRUE)

tax_columns <- colnames(mt_taxonomy_rar)[colnames(mt_taxonomy_rar) != "OTU_ID"]
mt_taxonomy_matrix <- mt_taxonomy_rar[, tax_columns]
mt_taxonomy_table_rar <- tax_table(as.matrix(mt_taxonomy_matrix))

# Match samples
matching_samples <- intersect(rownames(mt_metadata), 
                             colnames(mt_sbg_table_rar))
mt_metadata_filtered <- mt_metadata[matching_samples, , drop = FALSE]
mt_sbg_rar_filtered <- mt_sbg_table_rar[, matching_samples, drop = FALSE]

# Match OTUs
matching_otus <- intersect(rownames(mt_sbg_rar_filtered), 
                          rownames(mt_taxonomy_table_rar))
mt_sbg_rar_filtered <- mt_sbg_rar_filtered[matching_otus, ]
mt_taxonomy_table_rar_filtered <- mt_taxonomy_table_rar[matching_otus, ]

# Create phyloseq object
mt_physeq_rar <- phyloseq(
  otu_table(mt_sbg_rar_filtered, taxa_are_rows = TRUE),
  sample_data(mt_metadata_filtered),
  tax_table(mt_taxonomy_table_rar_filtered)
)

cat("\nPhyloseq object created successfully!\n")
cat("Number of taxa:", ntaxa(mt_physeq_rar), "\n")
cat("Number of samples:", nsamples(mt_physeq_rar), "\n")

# Save phyloseq object
saveRDS(mt_physeq_rar, "mt_physeq_rarefied.rds")

cat("\nTaxa abundance and qPCR analysis complete!\n")
cat("Results saved to working directory and 'results/' folder.\n")
