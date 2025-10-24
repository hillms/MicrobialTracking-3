################################################################################
# Antibiotic Resistance Genes (ARG), Virulence Factors (VF), 
# and Lateral Gene Transfer (LGT) Analysis
# MT-2 and MT-3 Microbiome Data
################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)
library(rstatix)
library(forcats)
library(hrbrthemes)
library(purrr)

# Define color palettes
colors_ARGs <- c("#009a9c", "#009a9c")
colors_ARGs_PMA <- c("#52B9BB", "#52B9BB", "#52B9BB")
colors_VFs <- c("#665AA6", "#665AA6")
colors_VFs_PMA <- c("#968EC2", "#968EC2")
colors_LGTs <- c('#4EA4D4', '#4EA4D4')
comparisons_Mission <- list(c("MT-2", "MT-3"))

################################################################################
# 1. DATA IMPORT AND PREPARATION
################################################################################

# Import metadata
MT <- read.delim(file = "MT-2_MT-3_combined_metadata_31Jan.txt", header = TRUE)
MT <- MT[MT$Flight != "8", ]

# Separate untreated and PMA-treated samples
MT_noPMA <- subset(MT, PMA_Treated == "No")
MT_PMA <- subset(MT, PMA_Treated == "Yes")

# Convert to numeric and remove NA surfaces
for(df_name in c("MT_noPMA", "MT_PMA")) {
  df <- get(df_name)
  df$LGTs_norm <- as.numeric(as.character(df$LGTs_norm))
  df$ARGs_norm <- as.numeric(as.character(df$ARGs_norm))
  df$VFs_norm <- as.numeric(as.character(df$VFs_norm))
  df <- df[!is.na(df$Surface), ]
  assign(df_name, df)
}

################################################################################
# 2. ARG ANALYSIS - UNTREATED SAMPLES
################################################################################

cat("\n=== ARG Analysis (Untreated) ===\n")

# Mission comparison
wilcox_results_nopma_mission_arg <- MT_noPMA %>%
  pairwise_wilcox_test(ARGs_norm ~ Mission, p.adjust.method = "BH")
print(wilcox_results_nopma_mission_arg)

# Flight comparison
wilcox_results_nopma_flight_arg <- MT_noPMA %>%
  pairwise_wilcox_test(ARGs_norm ~ Flight, p.adjust.method = "BH")
print(wilcox_results_nopma_flight_arg)

# Surface comparison (MT-3 only)
MT_noPMA_MT3 <- MT_noPMA[MT_noPMA$Mission != "MT-2", ]
wilcox_results_surface_nopma <- MT_noPMA_MT3 %>%
  pairwise_wilcox_test(ARGs_norm ~ Surface, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_surface_nopma)

# Collection day comparison
wilcox_results_day_nopma <- MT_noPMA_MT3 %>%
  pairwise_wilcox_test(ARGs_norm ~ Collection_Day, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_day_nopma)

################################################################################
# 3. ARG ANALYSIS - PMA-TREATED SAMPLES
################################################################################

cat("\n=== ARG Analysis (PMA-Treated) ===\n")

# Mission comparison
wilcox_results_pma_mission_arg <- MT_PMA %>%
  pairwise_wilcox_test(ARGs_norm ~ Mission, p.adjust.method = "BH")
print(wilcox_results_pma_mission_arg)

# Flight comparison
wilcox_results_pma_flight_arg <- MT_PMA %>%
  pairwise_wilcox_test(ARGs_norm ~ Flight, p.adjust.method = "BH")
print(wilcox_results_pma_flight_arg)

# Surface comparison (MT-3 only)
MT_PMA_MT3 <- MT_PMA[MT_PMA$Mission != "MT-2", ]
wilcox_results_surface_pma <- MT_PMA_MT3 %>%
  pairwise_wilcox_test(ARGs_norm ~ Surface, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_surface_pma)

################################################################################
# 4. VF ANALYSIS - UNTREATED SAMPLES
################################################################################

cat("\n=== VF Analysis (Untreated) ===\n")

# Mission comparison
wilcox_results_nopma_mission_vf <- MT_noPMA %>%
  pairwise_wilcox_test(VFs_norm ~ Mission, p.adjust.method = "BH")
print(wilcox_results_nopma_mission_vf)

# Flight comparison
wilcox_results_nopma_flight_vf <- MT_noPMA %>%
  pairwise_wilcox_test(VFs_norm ~ Flight, p.adjust.method = "BH")
print(wilcox_results_nopma_flight_vf)

# Surface comparison (MT-3 only)
wilcox_results_surface_nopma_vf <- MT_noPMA_MT3 %>%
  pairwise_wilcox_test(VFs_norm ~ Surface, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_surface_nopma_vf)

################################################################################
# 5. VF ANALYSIS - PMA-TREATED SAMPLES
################################################################################

cat("\n=== VF Analysis (PMA-Treated) ===\n")

# Mission comparison
wilcox_results_pma_mission_vf <- MT_PMA %>%
  pairwise_wilcox_test(VFs_norm ~ Mission, p.adjust.method = "BH")
print(wilcox_results_pma_mission_vf)

# Flight comparison (filter groups with variation)
MT_PMA_filtered <- MT_PMA %>%
  group_by(Flight) %>%
  filter(sum(VFs_norm > 0) > 1) %>%
  ungroup()

wilcox_results_pma_vf <- MT_PMA_filtered %>%
  pairwise_wilcox_test(VFs_norm ~ Flight, p.adjust.method = "BH")
print(wilcox_results_pma_vf)

################################################################################
# 6. LGT ANALYSIS - UNTREATED SAMPLES
################################################################################

cat("\n=== LGT Analysis (Untreated) ===\n")

# Flight comparison
wilcox_results_nopma_lgt <- MT_noPMA %>%
  pairwise_wilcox_test(LGTs_norm ~ Flight, p.adjust.method = "BH")
print(wilcox_results_nopma_lgt)

# Surface comparison (MT-3 only)
wilcox_results_surface_nopma_lgt <- MT_noPMA_MT3 %>%
  pairwise_wilcox_test(LGTs_norm ~ Surface, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_surface_nopma_lgt)

# Collection day comparison
wilcox_results_day_nopma_lgt <- MT_noPMA_MT3 %>%
  pairwise_wilcox_test(LGTs_norm ~ Collection_Day, p.adjust.method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"))
print(wilcox_results_day_nopma_lgt)

################################################################################
# 7. ARG MECHANISM ANALYSIS
################################################################################

cat("\n=== ARG Mechanism Analysis ===\n")

# Load mechanism data
arg_mech_mt2 <- read.delim("all_mech_MT2_MSH_13Nov2024.txt", header = TRUE)
arg_mech_mt3 <- read.delim("all_mech_MT3_MSH_13Nov2024.txt", header = TRUE)
arg_mech_metadata <- read.delim("MT-2_MT-3_combined_metadata_31Jan.txt", 
                                header = TRUE)

# Convert CPM to numeric and add mission labels
arg_mech_mt2$CPM <- as.numeric(arg_mech_mt2$CPM)
arg_mech_mt3$CPM <- as.numeric(arg_mech_mt3$CPM)
arg_mech_mt2$Mission <- "MT-2"
arg_mech_mt3$Mission <- "MT-3"

# Combine datasets
arg_mech_all <- bind_rows(arg_mech_mt2, arg_mech_mt3)

# Clean and standardize mechanism names
arg_mech_all <- arg_mech_all %>%
  mutate(
    Mechanism = ifelse(Mechanism %in% c("Multi-drug resistance", 
                                       "Multi-drug_resistance"),
                      "Multi-drug resistance", Mechanism),
    CPM = as.numeric(CPM)
  ) %>%
  group_by(Sample_ID, Mechanism, Mission) %>%
  summarise(CPM = sum(CPM, na.rm = TRUE), .groups = "drop")

# Merge with metadata
arg_merged_all <- arg_mech_all %>%
  left_join(arg_mech_metadata, by = "Sample_ID")

# Filter for untreated samples
arg_merged_noPMA <- arg_merged_all %>% 
  filter(PMA_Treated == "No")

# Calculate total CPM by mechanism
cat("\nTotal CPM by mechanism:\n")
print(aggregate(CPM ~ Mechanism, data = arg_merged_noPMA, FUN = sum, na.rm = TRUE))

# Statistical testing: Compare mechanisms between missions for each surface
wilcox_results_mech <- arg_merged_noPMA %>%
  group_by(Surface, Mechanism) %>%
  filter(n_distinct(Mission.x) == 2) %>%
  wilcox_test(CPM ~ Mission.x) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

wilcox_results_mech_sig <- wilcox_results_mech %>% 
  filter(p.adj <= 0.05)

cat("\nSignificant mechanism differences between missions:\n")
print(wilcox_results_mech_sig)

# Compare CPM between surfaces for each mechanism (MT-2)
mt2_surface_test <- arg_merged_noPMA %>%
  filter(Mission.x == "MT-2", !is.na(CPM)) %>%
  group_by(Mechanism) %>%
  filter(n_distinct(Surface) > 1, n() >= 4) %>%
  nest() %>%
  mutate(
    test_result = map(data, ~{
      tryCatch({
        wilcox.test(CPM ~ Surface, data = .x)
      }, error = function(e) NULL)
    })
  ) %>%
  filter(!map_lgl(test_result, is.null)) %>%
  mutate(p_value = map_dbl(test_result, ~.x$p.value)) %>%
  select(-data, -test_result) %>%
  mutate(
    p.adj = p.adjust(p_value, method = "BH"),
    Mission = "MT-2",
    p.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Compare CPM between surfaces for each mechanism (MT-3)
mt3_surface_test <- arg_merged_noPMA %>%
  filter(Mission.x == "MT-3", !is.na(CPM)) %>%
  group_by(Mechanism) %>%
  filter(n_distinct(Surface) > 1, n() >= 4) %>%
  nest() %>%
  mutate(
    test_result = map(data, ~{
      tryCatch({
        wilcox.test(CPM ~ Surface, data = .x)
      }, error = function(e) NULL)
    })
  ) %>%
  filter(!map_lgl(test_result, is.null)) %>%
  mutate(p_value = map_dbl(test_result, ~.x$p.value)) %>%
  select(-data, -test_result) %>%
  mutate(
    p.adj = p.adjust(p_value, method = "BH"),
    Mission = "MT-3",
    p.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("\nMT-2 - Surface differences by mechanism:\n")
print(mt2_surface_test)

cat("\nMT-3 - Surface differences by mechanism:\n")
print(mt3_surface_test)

# PMA-treated analysis
arg_merged_PMA <- arg_merged_all %>% 
  filter(PMA_Treated == "Yes")

wilcox_results_mech_pma <- arg_merged_PMA %>%
  group_by(Surface, Mechanism) %>%
  filter(n_distinct(Mission.x) == 2) %>%
  wilcox_test(CPM ~ Mission.x) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

wilcox_results_mech_sig_pma <- wilcox_results_mech_pma %>% 
  filter(p.adj <= 0.05)

cat("\nSignificant mechanism differences (PMA-treated):\n")
print(wilcox_results_mech_sig_pma)

################################################################################
# 8. EXPORT STATISTICAL RESULTS
################################################################################

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Export all statistical results
write.csv(wilcox_results_nopma_mission_arg, 
          "results/wilcox_ARG_mission_noPMA.csv", row.names = FALSE)
write.csv(wilcox_results_pma_mission_arg, 
          "results/wilcox_ARG_mission_PMA.csv", row.names = FALSE)
write.csv(wilcox_results_surface_nopma, 
          "results/wilcox_ARG_surface_noPMA.csv", row.names = FALSE)
write.csv(wilcox_results_surface_nopma_vf, 
          "results/wilcox_VF_surface_noPMA.csv", row.names = FALSE)
write.csv(wilcox_results_surface_nopma_lgt, 
          "results/wilcox_LGT_surface_noPMA.csv", row.names = FALSE)
write.csv(wilcox_results_mech_sig, 
          "results/wilcox_ARG_mechanism_significant_noPMA.csv", row.names = FALSE)
write.csv(wilcox_results_mech_sig_pma, 
          "results/wilcox_ARG_mechanism_significant_PMA.csv", row.names = FALSE)
write.csv(mt2_surface_test, 
          "results/mt2_mechanism_surface_comparison.csv", row.names = FALSE)
write.csv(mt3_surface_test, 
          "results/mt3_mechanism_surface_comparison.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved to 'results/' directory.\n")
