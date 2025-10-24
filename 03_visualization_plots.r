################################################################################
# Visualization and Plotting Scripts
# MT-2 and MT-3 Microbiome Data
################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)
library(hrbrthemes)
library(RColorBrewer)
library(pheatmap)

# Define color palettes
colors_ARGs <- c("#009a9c", "#009a9c")
colors_VFs <- c("#665AA6", "#665AA6")
colors_LGTs <- c('#4EA4D4', '#4EA4D4')

################################################################################
# 1. ARG/VF/LGT COMBINED PLOTS - PMA-TREATED SAMPLES
################################################################################

# Load data
MT <- read.delim(file = "MT-2_MT-3_combined_metadata_31Jan.txt", header = TRUE)
MT <- MT[MT$Flight != "8", ]
MT_PMA <- subset(MT, PMA_Treated == "Yes")
MT_PMA$LGTs_norm <- as.numeric(as.character(MT_PMA$LGTs_norm))
MT_PMA$ARGs_norm <- as.numeric(as.character(MT_PMA$ARGs_norm))
MT_PMA$VFs_norm <- as.numeric(as.character(MT_PMA$VFs_norm))
MT_PMA <- MT_PMA[!is.na(MT_PMA$Surface), ]

MT_noPMA <- subset(MT, PMA_Treated == "No")
MT_noPMA$LGTs_norm <- as.numeric(as.character(MT_noPMA$LGTs_norm))
MT_noPMA$ARGs_norm <- as.numeric(as.character(MT_noPMA$ARGs_norm))
MT_noPMA$VFs_norm <- as.numeric(as.character(MT_noPMA$VFs_norm))
MT_noPMA <- MT_noPMA[!is.na(MT_noPMA$Surface), ]

## ARG PLOTS ##

a1 <- ggplot(MT_PMA, aes(x = Mission, y = ARGs_norm, color = Mission, 
                         fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.5) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = colors_ARGs) +
  scale_fill_manual(values = colors_ARGs) +
  labs(x = "Mission", y = "AMR genes", tag = "A") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

b1 <- ggplot(MT_PMA, aes(x = factor(Flight), y = ARGs_norm, 
                         color = Mission, fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.5) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = colors_ARGs) +
  scale_fill_manual(values = colors_ARGs) +
  labs(x = "Flight", y = "AMR genes") +
  facet_wrap(~ Mission, scales = "free_x",
             labeller = as_labeller(c("MT-2" = "Mission 2",
                                     "MT-3" = "Mission 3"))) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

c1 <- ggplot(MT_PMA, aes(x = Mission, y = fct_rev(Surface), 
                         size = ARGs_norm, fill = ARGs_norm)) +
  geom_point(shape = 16, fill = "#009a9c", color = "#52B9BB") +
  scale_size(range = c(1, 12), name = "CPM") +
  labs(x = "Mission", y = "Surface") +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(hjust = 0.5),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

combined_plot_arg_PMA <- (a1 | b1 | c1) +
  plot_layout(widths = c(1.2, 2, 0.7)) +
  plot_annotation(
    title = "Antibiotic Resistance Gene (ARG) Proportion",
    theme = theme(plot.title = element_text(hjust = 0.5)))

## VF PLOTS ##

f <- ggplot(MT_PMA, aes(x = Mission, y = VFs_norm, color = Mission, 
                        fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = colors_VFs) +
  scale_fill_manual(values = colors_VFs) +
  labs(x = "Mission", y = "VFs", tag = "B") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

g <- ggplot(MT_PMA, aes(x = factor(Flight), y = VFs_norm, 
                        color = Mission, fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ Mission, scales = "free_x",
             labeller = as_labeller(c("MT-2" = "Mission 2",
                                     "MT-3" = "Mission 3"))) +
  scale_color_manual(values = colors_VFs) +
  scale_fill_manual(values = colors_VFs) +
  labs(x = "Flight", y = "VFs") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

h <- ggplot(MT_PMA, aes(x = Mission, y = fct_rev(Surface), 
                        size = VFs_norm, fill = VFs_norm)) +
  geom_point(shape = 16, fill = "#665AA6", color = "#665AA6") +
  scale_size(range = c(1, 12), name = "VF/Mbp") +
  labs(x = "Mission", y = "Surface") +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(hjust = 0.5),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

combined_plot_vf_PMA <- (f | g | h) +
  plot_layout(widths = c(1.2, 2, 0.7)) +
  plot_annotation(
    title = "Virulence Factor (VF) Proportion",
    theme = theme(plot.title = element_text(hjust = 0.5)))

## LGT PLOTS ##

lgt_1 <- ggplot(MT_noPMA, aes(x = Mission, y = LGTs_norm, 
                              color = Mission, fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(width = 0.2) +
  scale_color_manual(values = colors_LGTs) +
  scale_fill_manual(values = colors_LGTs) +
  labs(x = "Mission", y = "LGT events", tag = "C") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

lgt_2 <- ggplot(MT_noPMA, aes(x = factor(Flight), y = LGTs_norm, 
                              color = Mission, fill = Mission)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ Mission, scales = "free_x",
             labeller = as_labeller(c("MT-2" = "Mission 2",
                                     "MT-3" = "Mission 3"))) +
  scale_color_manual(values = colors_LGTs) +
  scale_fill_manual(values = colors_LGTs) +
  labs(x = "Flight", y = "LGT events") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

lgt_3 <- ggplot(MT_noPMA, aes(x = Mission, y = fct_rev(Surface), 
                              size = LGTs_norm, fill = LGTs_norm)) +
  geom_point(shape = 16, fill = "#4EA4D4", color = "#4EA4D4") +
  scale_size(range = c(1, 12), name = "LGTs/Mbp") +
  labs(x = "Mission", y = "Surface") +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(hjust = 0.5),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11)) +
  scale_x_discrete(labels = c("MT-2" = "Mission 2", "MT-3" = "Mission 3"))

combined_plot_lgt <- (lgt_1 | lgt_2 | lgt_3) +
  plot_layout(widths = c(1.2, 2, 0.7)) +
  plot_annotation(
    title = "Lateral Gene Transfer (LGT) events",
    theme = theme(plot.title = element_text(hjust = 0.5)))

## FINAL COMBINED PLOT ##

final_combined <- (combined_plot_arg_PMA / combined_plot_vf_PMA / combined_plot_lgt)

# Save plots
ggsave("figures/combined_ARG_VF_LGT.pdf", final_combined, 
       width = 14, height = 12)

################################################################################
# 2. ARG MECHANISM BUBBLE PLOTS
################################################################################

# Load mechanism data
arg_mech_mt2 <- read.delim("all_mech_MT2_MSH_13Nov2024.txt", header = TRUE)
arg_mech_mt3 <- read.delim("all_mech_MT3_MSH_13Nov2024.txt", header = TRUE)
arg_mech_metadata <- read.delim("MT-2_MT-3_combined_metadata_31Jan.txt", 
                                header = TRUE)

# Process MT-2 data
arg_mech_mt2 <- arg_mech_mt2 %>%
  mutate(
    Mechanism = ifelse(Mechanism %in% c("Multi-drug resistance", 
                                       "Multi-drug_resistance"),
                      "Multi-drug resistance", Mechanism),
    CPM = as.numeric(CPM)) %>%
  group_by(Sample_ID, Mechanism) %>%
  summarise(CPM = sum(CPM), .groups = "drop")

arg_merged_mt2 <- arg_mech_mt2 %>%
  left_join(arg_mech_metadata, by = "Sample_ID")

arg_merged_mt2_noPMA <- subset(arg_merged_mt2, PMA_Treated == "No")

arg_merged_mt2_agg <- arg_merged_mt2_noPMA %>%
  group_by(Surface, Mechanism) %>%
  summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

e.a <- ggplot(arg_merged_mt2_agg, aes(x = Surface, y = Mechanism, size = CPM)) +
  geom_point(aes(size = CPM), color = "#009a9c", fill = "#009a9c", 
             shape = 16, show.legend = TRUE) +
  scale_size(range = c(1, 11), limits = c(0, 4000)) +
  theme_minimal() +
  labs(x = "Surface", y = "Mechanism", size = "CPM", title = "MT-2", tag = "A") +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11))

# Process MT-3 data
arg_mech_mt3 <- arg_mech_mt3 %>%
  mutate(
    Mechanism = ifelse(Mechanism %in% c("Multi-drug resistance", 
                                       "Multi-drug_resistance"),
                      "Multi-drug resistance", Mechanism),
    CPM = as.numeric(CPM)) %>%
  group_by(Sample_ID, Mechanism) %>%
  summarise(CPM = sum(CPM), .groups = "drop")

arg_merged_mt3 <- arg_mech_mt3 %>%
  left_join(arg_mech_metadata, by = "Sample_ID")

arg_merged_mt3_noPMA <- subset(arg_merged_mt3, PMA_Treated == "No")

arg_merged_mt3_agg <- arg_merged_mt3_noPMA %>%
  group_by(Surface, Mechanism) %>%
  summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

e.b <- ggplot(arg_merged_mt3_agg, aes(x = Surface, y = Mechanism, size = CPM)) +
  geom_point(aes(size = CPM), color = "#009a9c", fill = "#009a9c", 
             shape = 16, show.legend = TRUE) +
  scale_size(range = c(1, 11), limits = c(0, 4000)) +
  theme_minimal() +
  labs(x = "Surface", y = "Mechanism", size = "CPM", title = "MT-3") +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11))

arg_mech_plot_combine <- (e.a | e.b)

ggsave("figures/ARG_mechanism_bubbleplot.pdf", arg_mech_plot_combine, 
       width = 12, height = 6)

################################################################################
# 3. ARG RESISTANCE HEATMAP
################################################################################

# Read data
arg_heatmap <- read.delim(file = "Isolate_AMR_Heatmap.txt", header = TRUE)
arg_heatmap <- arg_heatmap %>% arrange(Species)
rownames(arg_heatmap) <- arg_heatmap$Strain_ID

# Define antibiotic classes
antibiotic_classes <- data.frame(
  Antibiotic = c("Amoxicillin_ClavulanicAcid", "Aztreonam", "Cefoxitin", 
                 "Ceftriaxone", "Cefazolin", "Cloxacillin", "Meropenem", 
                 "Oxacillin", "Penicillin", "Ciprofloxacin", "Gemifloxacin", 
                 "Levofloxacin", "Moxifloxacin", "Norfloxacin", "Ofloxacin",
                 "Clindamycin", "Azithromycin", "Erythromycin"),
  Class = c(rep("Beta-lactam", 9), rep("Fluoroquinolone", 6),
            "Lacosamide", rep("Macrolide", 2))
)

# Extract available antibiotics
available_antibiotics <- intersect(antibiotic_classes$Antibiotic, 
                                  colnames(arg_heatmap))
antibiotic_classes <- antibiotic_classes[
  antibiotic_classes$Antibiotic %in% available_antibiotics, ]

# Order antibiotics by class
ordered_abx <- antibiotic_classes %>%
  arrange(Class, Antibiotic) %>%
  pull(Antibiotic)

abx_matrix <- arg_heatmap[, ordered_abx]

# Create annotations
col_ann <- data.frame(Class = antibiotic_classes$Class)
rownames(col_ann) <- antibiotic_classes$Antibiotic

# Define colors
resistance_colors <- colorRampPalette(c("white", "#F6C87D"))(2)
custom_colors_anti_class <- c(
  "Beta-lactam" = "#525174",
  "Fluoroquinolone" = "dodgerblue2",
  "Lacosamide" = "#5DD39E",
  "Macrolide" = "#BCE784"
)

ann_colors <- list(Class = custom_colors_anti_class)

# Create italicized species names
species_names <- arg_heatmap$Species
italicized_species <- lapply(species_names, function(x) {
  parts <- strsplit(x, " ")[[1]]
  if(length(parts) >= 2) {
    bquote(italic(.(parts[1])) ~ italic(.(parts[2])))
  } else {
    bquote(italic(.(x)))
  }
})

italicized_species_expr <- do.call(expression, italicized_species)
rownames(abx_matrix) <- species_names

# Create heatmap
pdf("figures/ARG_resistance_heatmap.pdf", width = 10, height = 8)
pheatmap(
  mat = abx_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = col_ann,
  annotation_colors = ann_colors,
  color = resistance_colors,
  show_rownames = TRUE,
  fontsize_row = 10,
  labels_row = italicized_species_expr,
  gaps_col = cumsum(table(factor(antibiotic_classes$Class, 
                                levels = unique(antibiotic_classes$Class))))[
    1:(length(unique(antibiotic_classes$Class))-1)],
  legend = FALSE,
  angle_col = 315
)
dev.off()

################################################################################
# 4. VIRULENCE SURVIVAL CURVES
################################################################################

# Load virulence data
VF_all <- read.delim(file = "most_virulent_all_missions.txt", header = TRUE)

custom_colors_vf <- c("#f46d43", "#fdae61", "#fee090", "#74add1", "#4A6A80")
custom_order <- c("Escherichia coli", "Bacillus cereus", 
                 "Pseudomonas granadensis", "Staphylococcus saprophyticus",
                 "Pseudomonas aeruginosa")

VF_all$Taxon <- factor(VF_all$Taxon, levels = custom_order)

vf_plot <- ggplot(VF_all, aes(x = Hours, y = Survival, color = Taxon, 
                              shape = Treatment)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 24, 30, 48)) +
  labs(x = "Hours", y = "Percent survival", color = "Taxon", 
       shape = "Treatment") +
  scale_color_manual(values = custom_colors_vf) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))

ggsave("figures/virulence_survival_curves.pdf", vf_plot, 
       width = 10, height = 6)

################################################################################
# 5. PREDICTED PHENOTYPE BAR PLOTS
################################################################################

pheno <- read.delim(file = "percent_pheno.txt", header = TRUE)

# Normalize within each Assay
df_norm <- pheno %>%
  group_by(Assay) %>%
  mutate(Percent = Abundance / sum(Abundance) * 100)

# Define colors
phenotype_colors <- c(
  "FP" = "#66c2a5",
  "FN" = "#F6C87D",
  "Predicted" = "#8da0cb",
  "Both" = "#e78ac3"
)

# AMR plot
plot_AMR <- ggplot(filter(df_norm, Assay == "AMR"),
                   aes(x = Assay, y = Percent, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 4) +
  labs(title = "AMR", y = "Percent", x = NULL, tag = "A") +
  theme_minimal() +
  scale_fill_manual(values = phenotype_colors, drop = FALSE)

# Virulence plot
plot_Virulence <- ggplot(filter(df_norm, Assay == "Virulence"),
                         aes(x = Assay, y = Percent, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 4) +
  labs(title = "Virulence", y = "Percent", x = NULL, tag = "B") +
  theme_minimal() +
  scale_fill_manual(values = phenotype_colors, drop = FALSE)

pheno_combined <- (plot_AMR | plot_Virulence)
ggsave("figures/predicted_phenotype.pdf", pheno_combined, 
       width = 8, height = 5)

cat("\nAll plots generated and saved to 'figures/' directory!\n")