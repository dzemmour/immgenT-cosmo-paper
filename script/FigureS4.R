# Extended Data Figure 4. Cluster robustness metrics.
#
# Panels produced:
#   S4a  Experiment composition of each cluster.
#   S4b  Cluster proportions across the spleen standard samples.
#   S4d  Cluster separability score per cluster across experiments, in protein
#        space.
#
# Panel S4c (cross-experiment cosine similarity of each cluster's expression
# profile) is not reproduced here: the similarity matrix is not among the cached
# inputs. See analysis/FigureS4.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "S4: Cluster robustness" (line 1498),
# subsections "S4a: Clusters by IGT" and "S4d: Cluster Separability in ADT
# Space", plus the "Boxplot_Level2ByIGT_spleencontrols" chunk (line 1976), which
# is panel b and sits further down the file under the composition analysis.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                         [primary input]
#   color_palette_igt.csv                                     [curated input]
#   ClusterSeparation_ADT_IGT.Rds                             [code/pipeline/01]
#   annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds [primary input]
#   Sample_metadata_david_20260107_v11.csv                    [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
})

source("code/R/setup.R")
source("code/R/composition.R")

figure_dir <- "Extended Data Figure 4"

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus") %>% rownames()]
mypal_igt <- load_igt_palette()

# ============================================================
# S4a: Experiment composition of each cluster
# ============================================================
# Bars sum to 1 within a cluster, after each experiment's contribution has
# first been normalised by its own size (see CompositionOfEachCluster). A
# cluster made of one colour would be experiment-specific and therefore
# suspect; the panel's claim is that none is.
comp_df <- CompositionOfEachCluster(so@meta.data,
                                    cluster_col = "annotation_level2",
                                    sample_col = "IGT")

comp_df$annotation_level1 <- so_orig$annotation_level1[
    match(comp_df$annotation_level2, so_orig$annotation_level2)
]
comp_df$annotation_level1 <- factor(comp_df$annotation_level1,
                                    levels = levels(so$annotation_level1))

p_s4a <- ggplot(comp_df) +
    geom_bar(aes(x = label, y = p_weighted, fill = IGT), stat = "identity") +
    scale_fill_manual(values = mypal_igt) +
    facet_grid(. ~ annotation_level1, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Experiment composition") +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold")
    )

panel_pdf(figure_dir, "S4a_cluster_composition_by_IGT", 15, 5)
print(p_s4a + NoGrid() + NoLegend())
dev.off()

# ============================================================
# S4b: Cluster proportions across the spleen standard samples
# ============================================================
# The technical replicate: the same standardised spleen preparation, run in
# many experiments. Because the input material is nominally identical, the
# spread of each cluster's proportion is a direct read-out of measurement
# reproducibility rather than of biology.
mdata <- load_sample_composition(data_path)

# Two spleen samples that belong to the standard series but were not flagged as
# such in the metadata export.
so_orig$spleen_standard[
    so_orig$sample_code %in% c("I44H2_spleen_allT_mouse_0162",
                               "I45H1_spleen_CD45p_mouse_0170")
] <- TRUE
mdata$spleen_standard <-
    mdata$IGTHT %in% unique(so_orig$IGTHT[so_orig$spleen_standard])

df_long <- sample_composition_long(
    mdata %>% filter(spleen_standard),
    levels(so_orig$annotation_level2),
    keep_ncells = TRUE
) %>%
    # Clusters represented by fewer than 10 cells in a sample give a proportion
    # too noisy to compare across replicates.
    filter(ncells > 10)

df_long$IGT <- so_orig$IGT[match(df_long$IGTHT, so_orig$IGTHT)]

p_s4b <- ggplot(df_long,
                aes(x = annotation_level2, y = prop,
                    fill = annotation_level1)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.5) +
    geom_jitter(aes(color = IGT), width = 0.15, height = 0,
                size = 1.2, alpha = 0.8) +
    scale_fill_manual(values = mypal_level1) +
    scale_color_manual(values = mypal_igt) +
    facet_grid(. ~ annotation_level1, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Proportion across experiments") +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold")
    ) +
    NoGrid() + NoLegend()

panel_pdf(figure_dir, "S4b_cluster_proportions_spleen_standards", 15, 5)
print(p_s4b)
dev.off()

# ============================================================
# S4d: Cluster separability in protein space
# ============================================================
# The counterpart of Figure 4c, computed on principal components of the
# CITE-seq protein matrix instead of the transcriptome. Clusters were defined
# on RNA, so a positive score here means the RNA-defined partition is also
# visible in an independent measurement modality.
sep_list <- readRDS(sprintf("%s/ClusterSeparation_ADT_IGT.Rds", data_path))
sep <- Reduce(rbind, sep_list)

dict_l1 <- so@meta.data %>% count(annotation_level1, annotation_level2)
sep$annotation_level1 <- dict_l1$annotation_level1[
    match(sep$annotation_level2, dict_l1$annotation_level2)
]
sep$annotation_level2 <- factor(sep$annotation_level2,
                                levels = levels(so$annotation_level2))

p_s4d <- ggplot(sep, aes(x = annotation_level2, y = s)) +
    geom_jitter(aes(color = IGT), width = 0.15, height = 0,
                size = 1.8, alpha = 0.9) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "black",
                 width = 0.6) +
    scale_color_manual(values = mypal_igt) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.3) +
    facet_grid(. ~ annotation_level1, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Cluster separability score") +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold")
    )

panel_pdf(figure_dir, "S4d_cluster_separability_ADT", 25, 6)
print(p_s4d + NoLegend())
dev.off()
