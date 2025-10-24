################################################################################
# Differential Abundance Analysis using ALDEx2
# MT-2 and MT-3 Microbiome Data
################################################################################

# Load required libraries
library(ALDEx2)
library(phyloseq)
library(dplyr)
library(ggplot2)

################################################################################
# 1. MISSION COMPARISON (Untreated Samples)
################################################################################

# Load data
feature_table <- read.table("physeq_feature_table_mt_unrar.txt", 
                            header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("MT-2_MT-3_combined_metadata_31Jan.txt", 
                       header = TRUE, row.names = 1, sep = "\t")
mt_taxonomy_unrar <- read.delim("physeq_taxonomy_mt_unrar.txt", 
                                header = TRUE, sep = "\t")
row.names(mt_taxonomy_unrar) <- mt_taxonomy_unrar$OTU_ID

# Filter for untreated samples only
metadata_nopma <- metadata[metadata$PMA_Treated != "Yes", ]

# Identify common samples
common_samples <- intersect(colnames(feature_table), rownames(metadata_nopma))
feature_table_nopma <- feature_table[, common_samples]
metadata_nopma <- metadata_nopma[common_samples, ]

# Check alignment
stopifnot(all(colnames(feature_table_nopma) == rownames(metadata_nopma)))

# Define groups
group <- metadata_nopma$Mission

# Run ALDEx2
aldex_results_nopma <- aldex(feature_table_nopma, 
                             conditions = group, 
                             mc.samples = 128, 
                             test = "t", 
                             effect = TRUE, 
                             denom = "iqlr")

# Add taxonomy
aldex_results_with_taxa <- merge(aldex_results_nopma, mt_taxonomy_unrar, 
                                 by = "row.names", all.x = TRUE)
colnames(aldex_results_with_taxa)[1] <- "OTU_ID"

# Filter significant taxa
significant_taxa <- aldex_results_with_taxa[aldex_results_with_taxa$we.ep < 0.05, ]

# Export results
write.csv(significant_taxa, "significant_taxa_MT2_MT3_nopma.csv", row.names = FALSE)
write.table(aldex_results_nopma, "aldex2_results_nopma.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)

# Visualize
plot_data <- significant_taxa %>%
  dplyr::select(Species, effect) %>%
  mutate(abs_effect = abs(effect)) %>%
  arrange(effect) %>%
  head(25) %>%
  mutate(Species = factor(Species, levels = unique(Species)))

ggplot(plot_data, aes(y = Species, x = "MT-3", fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = round(effect, 2)), color = "black") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Taxon", fill = "Effect")

################################################################################
# 2. MISSION COMPARISON (PMA-Treated Samples)
################################################################################

metadata_mission_pma <- metadata[metadata$PMA_Treated != "No", ]
common_samples_mission_pma <- intersect(colnames(feature_table), 
                                        rownames(metadata_mission_pma))
feature_table_mission_pma <- feature_table[, common_samples_mission_pma]
metadata_mission_pma <- metadata_mission_pma[common_samples_mission_pma, ]

stopifnot(all(colnames(feature_table_mission_pma) == rownames(metadata_mission_pma)))

group_mission <- metadata_mission_pma$Mission

aldex_results_mission_pma <- aldex(feature_table_mission_pma, 
                                   conditions = group_mission, 
                                   mc.samples = 128, 
                                   test = "t", 
                                   effect = TRUE, 
                                   denom = "iqlr")

aldex_results_mission_pma_with_taxa <- merge(aldex_results_mission_pma, 
                                             mt_taxonomy_unrar, 
                                             by = "row.names", all.x = TRUE)
colnames(aldex_results_mission_pma_with_taxa)[1] <- "OTU_ID"

significant_taxa_mission_pma <- aldex_results_mission_pma_with_taxa[
  aldex_results_mission_pma_with_taxa$we.ep < 0.05, ]

write.csv(significant_taxa_mission_pma, "significant_taxa_MT2_MT3_pma.csv", 
          row.names = FALSE)

################################################################################
# 3. PMA vs UNTREATED COMPARISON (MT-3 only)
################################################################################

metadata_pma_untreated <- metadata[metadata$Mission != "MT-2", ]
common_samples_pma_untreated <- intersect(colnames(feature_table), 
                                          rownames(metadata_pma_untreated))
feature_table_pma_untreated <- feature_table[, common_samples_pma_untreated]
metadata_pma_untreated <- metadata_pma_untreated[common_samples_pma_untreated, ]

stopifnot(all(colnames(feature_table_pma_untreated) == 
              rownames(metadata_pma_untreated)))

group_pma <- factor(metadata_pma_untreated$PMA_Treated, levels = c("Yes", "No"))

aldex_results_pma_untreated <- aldex(feature_table_pma_untreated, 
                                     conditions = group_pma, 
                                     mc.samples = 128, 
                                     test = "t", 
                                     effect = TRUE, 
                                     denom = "iqlr")

aldex_results_pma_untreated_with_taxa <- merge(aldex_results_pma_untreated, 
                                               mt_taxonomy_unrar, 
                                               by = "row.names", all.x = TRUE)
colnames(aldex_results_pma_untreated_with_taxa)[1] <- "OTU_ID"

significant_taxa_pma_untreated <- aldex_results_pma_untreated_with_taxa[
  aldex_results_pma_untreated_with_taxa$we.ep < 0.05, ]

write.csv(significant_taxa_pma_untreated, "significant_taxa_pma_untreated.csv", 
          row.names = FALSE)

################################################################################
# 4. SURFACE COMPARISON (MT-3, PMA-treated, pairwise vs Cupola)
################################################################################

metadata_surface_nopma <- metadata[metadata$Mission != "MT-2", ]
metadata_surface_nopma <- metadata_surface_nopma[
  metadata_surface_nopma$PMA_Treated != "No", ]
metadata_surface_nopma <- metadata_surface_nopma[
  metadata_surface_nopma$Surface != "NA", ]

common_samples_surface_nopma <- intersect(colnames(feature_table), 
                                          rownames(metadata_surface_nopma))
feature_table_surface_nopma <- feature_table[, common_samples_surface_nopma]
metadata_surface_nopma <- metadata_surface_nopma[common_samples_surface_nopma, ]

stopifnot(all(colnames(feature_table_surface_nopma) == 
              rownames(metadata_surface_nopma)))

group_surface_nopma <- factor(metadata_surface_nopma$Surface, 
                              levels = c("Cupola", "ARED", "Crew quarters", 
                                       "Dining table", "LAB", "PMM", 
                                       "Stowage rack", "WHC"))

# Pairwise comparisons against Cupola
results_list <- list()
other_surfaces <- setdiff(levels(group_surface_nopma), "Cupola")

for(surface in other_surfaces) {
  current_samples <- rownames(metadata_surface_nopma)[
    metadata_surface_nopma$Surface %in% c("Cupola", surface)]
  current_features <- feature_table[, current_samples]
  current_metadata <- metadata_surface_nopma[current_samples, ]
  
  current_group <- factor(current_metadata$Surface, 
                         levels = c("Cupola", surface))
  
  results_list[[surface]] <- aldex(current_features,
                                  conditions = current_group,
                                  mc.samples = 128,
                                  test = "t",
                                  effect = TRUE,
                                  denom = "all")
}

# Combine results
all_results <- data.frame(OTU_ID = character(), 
                          Surface = character(), 
                          effect = numeric(), 
                          p_value = numeric(), 
                          stringsAsFactors = FALSE)

for (surface in names(results_list)) {
  current_result <- results_list[[surface]]
  temp_df <- data.frame(
    OTU_ID = rownames(current_result),
    Surface = surface,
    effect = current_result$effect,
    p_value = current_result$we.ep,
    stringsAsFactors = FALSE
  )
  all_results <- rbind(all_results, temp_df)
}

# Add taxonomy
all_results_with_taxa <- merge(all_results, mt_taxonomy_unrar, 
                               by.x = "OTU_ID", by.y = 0, all.x = TRUE)
significant_results <- all_results_with_taxa[all_results_with_taxa$p_value < 0.05, ]

# Get top 50 species
species_max_effect <- significant_results %>%
  group_by(Species) %>%
  summarize(max_abs_effect = max(abs(effect), na.rm = TRUE)) %>%
  arrange(desc(max_abs_effect)) %>%
  head(50)

heatmap_data <- significant_results %>%
  filter(Species %in% species_max_effect$Species)

# Create heatmap
ggplot(heatmap_data, aes(x = Surface, y = Species, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = round(effect, 2)), color = "black", size = 3) +
  geom_hline(yintercept = seq(1.5, length(unique(heatmap_data$Species))), 
             color = "gray90", linewidth = 0.2) +
  geom_vline(xintercept = seq(1.5, length(unique(heatmap_data$Surface))), 
             color = "gray90", linewidth = 0.2) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_discrete(limits = rev) +
  labs(y = "Taxon", x = "Surface (compared to Cupola)", fill = "Effect Size")

################################################################################
# 5. FLIGHT COMPARISON (MT-3, PMA-treated)
################################################################################

metadata_day_nopma <- metadata[metadata$Mission != "MT-2", ]
metadata_day_nopma <- metadata_day_nopma[metadata_day_nopma$PMA_Treated != "No", ]
metadata_day_nopma <- metadata_day_nopma[metadata_day_nopma$Surface != "NA", ]

common_samples_day_nopma <- intersect(colnames(feature_table), 
                                      rownames(metadata_day_nopma))
feature_table_day_nopma <- feature_table[, common_samples_day_nopma]
metadata_day_nopma <- metadata_day_nopma[common_samples_day_nopma, ]

stopifnot(all(colnames(feature_table_day_nopma) == rownames(metadata_day_nopma)))

group_day_nopma <- factor(metadata_day_nopma$Flight, levels = c("9", "10", "11"))

# Pairwise comparisons against Flight 9
results_list_flight <- list()
other_days <- setdiff(levels(group_day_nopma), "9")

for(flight in other_days) {
  current_samples_day <- rownames(metadata_day_nopma)[
    metadata_day_nopma$Flight %in% c("9", flight)]
  current_features_day <- feature_table_day_nopma[, current_samples_day]
  current_metadata_day <- metadata_day_nopma[current_samples_day, ]
  
  current_group_day <- factor(current_metadata_day$Flight, 
                              levels = c("9", flight))
  
  results_list_flight[[flight]] <- aldex(current_features_day,
                                        conditions = current_group_day,
                                        mc.samples = 128,
                                        test = "t",
                                        effect = TRUE,
                                        denom = "all")
}

# Process results
aldex_df_list <- list()
for(flight in other_days) {
  aldex_result <- results_list_flight[[flight]]
  aldex_df <- data.frame(aldex_result)
  aldex_df$OTU_ID <- rownames(aldex_df)
  aldex_df$Flight <- flight
  aldex_df <- aldex_df[, c("OTU_ID", "Flight", "effect")]
  aldex_df_list[[flight]] <- aldex_df
}

heatmap_data_day_nopma <- do.call(rbind, aldex_df_list)
aldex_results_day_nopma_with_taxa <- merge(heatmap_data_day_nopma, 
                                           mt_taxonomy_unrar, 
                                           by = "OTU_ID", all.x = TRUE)

significant_taxa_day_nopma <- unique(
  aldex_results_day_nopma_with_taxa$Species[
    abs(aldex_results_day_nopma_with_taxa$effect) > 0.5])

heatmap_data_significant_day_nopma <- aldex_results_day_nopma_with_taxa[
  aldex_results_day_nopma_with_taxa$Species %in% significant_taxa_day_nopma, ]

write.table(heatmap_data_significant_day_nopma,
            file = "significant_taxa_MT3_flight_pma.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nDifferential abundance analysis complete!\n")
cat("Results saved to working directory.\n")
