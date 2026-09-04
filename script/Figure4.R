# Figure 4. 107 immgenT consensus clusters preserve heterogeneity and enhance
# biological interpretability.
#
# Panels produced:
#   4a  Saturation: cumulative number of clusters reaching 50 cells as
#       experiments accumulate.
#   4b  Cluster silhouette width before vs. after integration.
#   4c  Cluster separability score per cluster across experiments (RNA space).
#   4d  One sample's cluster proportions against the distribution across all
#       samples.
#   4e  That sample's CD8 cells on the CD8-specific MDE.
#   4f  That sample's experiment-specific UMAP, CD8.Q cells highlighted.
#   4g  That sample's CD8.Q cells in CD103 vs ITGB7 protein space.
#   4h-4j  The same three views for a naive small-intestine sample.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Figure 4a: cumulative number of
# clusters" (line 1274), "Fig 4b- Biological conservation analysis" (line 1346),
# "Fig 4c - Cluster Separability in RNA Space" (line 1395), "Composition
# analysis Figure 4" (line 1760), "Specific IGT examples with CD8 cl10/CD8.Q"
# (line 1885) and "ImmgenT clusters in IGT UMAPs" (line 1892).
#
# Panels 4e/4h come from the chunk DZ added at cosmo_paper.Rmd line 1946, which
# is a sketch rather than working code: it passes lineage = "CD4" (this figure
# is about CD8), references an undefined `bg_cl` layer left over from the
# ablation section, assigns `s` twice in a row so only the second sample
# survives, and builds `p2` without printing it. Reconstructed here as the two
# panels the published figure shows, restricted to the sample's CD8 cells so
# that the cells drawn are on the same lineage embedding as the background.
#
# 4b and 4c read cached per-experiment metrics rather than recomputing them;
# each needs an independent clustering or PCA of every experiment. See
# code/pipeline/01_cluster_metrics.R and code/R/separability.R.
#
# The Rmd also plotted the I21H9 composition boxplot (the counterpart of 4d for
# the naive sample), which is not published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                        [primary input]
#   color_palette_igt.csv                                    [curated input]
#   BioConservationSilhouette_IGT.Rds                        [code/pipeline/01]
#   ClusterSeparation_RNA_IGT.Rds                            [code/pipeline/01]
#   annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds [primary input]
#   Sample_metadata_david_20260107_v11.csv                   [curated input]
#   IGT40_seurat.Rds, IGT21_seurat.Rds                       [primary inputs]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source("code/R/setup.R")
source("code/R/utils.R")
source("code/R/composition.R")

figure_dir <- "Figure 4"
set.seed(1)

# The two samples shown in panels d-j.
SAMPLE_LCMV  <- "I40H8_smallintestineLP_LCMVarm_D60_allT_P14_SMARTA_mouse0142"
SAMPLE_NAIVE <- "I21H9_smallintestineLP_allT_mouse0036"

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus", cite_seq == TRUE) %>%
                  rownames()]
mypal_igt <- load_igt_palette()

# ============================================================
# 4a: Saturation of cluster discovery across experiments
# ============================================================
# A cluster counts as discovered at the experiment where its running cell total
# first exceeds 50. Without a threshold a single cell assigned in the first
# experiment would count as a discovery and the curve would saturate at once.
df_counts <- so@meta.data %>%
    as.data.frame() %>%
    count(IGT, annotation_level2, name = "n_cells") %>%
    mutate(IGT_num = as.numeric(gsub("[^0-9]", "", IGT)))

df_cum <- df_counts %>%
    arrange(annotation_level2, IGT_num) %>%
    group_by(annotation_level2) %>%
    mutate(cum_n_cells = cumsum(n_cells)) %>%
    ungroup()

discovery_tbl <- df_cum %>%
    filter(cum_n_cells > 50) %>%
    group_by(annotation_level2) %>%
    summarise(first_IGT_num = min(IGT_num), .groups = "drop")

igt_steps <- df_cum %>%
    distinct(IGT, IGT_num) %>%
    arrange(IGT_num) %>%
    rowwise() %>%
    mutate(n_new = sum(discovery_tbl$first_IGT_num == IGT_num)) %>%
    ungroup() %>%
    mutate(cum_clusters = cumsum(n_new))

igt_steps$IGT <- factor(igt_steps$IGT, levels = igt_steps$IGT)

p_4a <- ggplot(igt_steps, aes(x = IGT, y = cum_clusters, group = 1)) +
    geom_line() +
    geom_point(size = 1.5) +
    # 107 is the final cluster count.
    geom_hline(yintercept = 107, colour = "brown", linetype = "dashed",
               linewidth = 0.5) +
    labs(x = "IGT", y = "Cumulative # clusters (> 50 cells)") +
    ylim(0, 115) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

panel_pdf(figure_dir, "4a_saturation_clusters", 10, 5)
print(p_4a + NoGrid())
dev.off()

# ============================================================
# 4b: Cluster silhouette width, before vs. after integration
# ============================================================
sil_summary_list <- readRDS(sprintf("%s/BioConservationSilhouette_IGT.Rds",
                                    data_path))
sil_summary_merged <- Reduce(rbind, sil_summary_list)

fit <- lm(mean_sil_width.post ~ mean_sil_width.pre, data = sil_summary_merged)
eq_lab <- sprintf("y = %.2f + %.2f·x\nR² = %.2f",
                  coef(fit)[1], coef(fit)[2], summary(fit)$r.squared)

p_4b <- ggplot(sil_summary_merged) +
    geom_point(aes(mean_sil_width.pre, mean_sil_width.post, color = IGT),
               size = 1) +
    geom_smooth(aes(mean_sil_width.pre, mean_sil_width.post),
                method = "lm", se = FALSE, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black",
               linewidth = 0.3) +
    scale_color_manual(values = mypal_igt) +
    xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
    labs(x = "Mean silhouette width pre-integration",
         y = "Mean silhouette width post-integration") +
    theme_bw() + NoGrid() +
    annotate("text",
             x = min(sil_summary_merged$mean_sil_width.pre, na.rm = TRUE),
             y = max(sil_summary_merged$mean_sil_width.post, na.rm = TRUE),
             hjust = 0, vjust = 1, label = eq_lab, size = 3)

panel_pdf(figure_dir, "4b_silhouette_pre_post_integration", 5, 5)
print(p_4b + NoLegend())
dev.off()

# ============================================================
# 4c: Cluster separability per cluster (RNA space)
# ============================================================
sep_list <- readRDS(sprintf("%s/ClusterSeparation_RNA_IGT.Rds", data_path))
sep <- Reduce(rbind, sep_list)

dict_l1 <- so@meta.data %>% count(annotation_level1, annotation_level2)
sep$annotation_level1 <- dict_l1$annotation_level1[
    match(sep$annotation_level2, dict_l1$annotation_level2)
]
sep$annotation_level2 <- factor(sep$annotation_level2,
                                levels = levels(so$annotation_level2))

p_4c <- ggplot(sep, aes(x = annotation_level2, y = s)) +
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

panel_pdf(figure_dir, "4c_cluster_separability_RNA", 25, 6)
print(p_4c + NoLegend())
dev.off()

# ============================================================
# 4d: One sample's cluster proportions against all samples
# ============================================================
# Boxplots are the distribution across every sample in the atlas; the coloured
# points are this one sample.
mdata <- load_sample_composition(data_path)
df_long <- sample_composition_long(mdata, levels(so_orig$annotation_level2))

p_boxes <- ggplot(df_long, aes(x = variable, y = value)) +
    geom_boxplot(outlier.shape = NA, fill = "grey", alpha = 0.5)

example <- df_long %>% filter(IGTHT == "I40H8")

p_4d <- p_boxes +
    geom_point(data = example,
               aes(variable, value, color = variable), size = 2, alpha = 1) +
    facet_grid(. ~ annotation_level1, scales = "free_x", space = "free_x") +
    scale_color_manual(values = mypal_level2) +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) +
    labs(x = NULL, y = "Proportion (of each lineage)") +
    ggtitle(SAMPLE_LCMV)

panel_pdf(figure_dir, "4d_sample_composition_I40H8", 15, 5)
print(p_4d + NoLegend() + NoGrid())
dev.off()

# ============================================================
# 4e/4h: Sample CD8 cells on the CD8-specific MDE
# ============================================================
# All immgenT CD8 cells form the background, as a pre-rendered raster rather
# than ~200,000 plotted points, which keeps the PDF openable. The sample's own
# CD8 cells are drawn over it as labelled points.
bg_cd8 <- AddEmbeddingRasterBackground(
    lineage = "CD8", interpolate = FALSE, color = "black", alpha = 0.5
)

cd8_mde_panel <- function(sample_code, title) {
    cells <- so_orig@meta.data %>%
        filter(sample_code == !!sample_code, annotation_level1 == "CD8") %>%
        rownames()
    p <- DimPlot(
        so_orig[, cells],
        reduction = "mde_incremental",
        group.by = "annotation_level2",
        pt.size = 2, label = TRUE, raster = TRUE,
        cols = mypal_level2
    )
    p$layers <- c(list(bg_cd8), p$layers)
    p +
        theme_void() +
        coord_fixed(xlim = c(-2.5, 3.5), ylim = c(-2.5, 2.5)) +
        ggtitle(title) +
        NoLegend()
}

panel_pdf(figure_dir, "4e_CD8MDE_I40H8", 6, 6)
print(cd8_mde_panel(
    SAMPLE_LCMV,
    "CD8+ T cells, small intestine lamina propria, 30 days post LCMV arm"
))
dev.off()

panel_pdf(figure_dir, "4h_CD8MDE_I21H9", 6, 6)
print(cd8_mde_panel(
    SAMPLE_NAIVE,
    "CD8+ T cells, small intestine lamina propria, naive mouse"
))
dev.off()

# ============================================================
# 4f/4i: Experiment-specific UMAP with CD8.Q highlighted
# ============================================================
# Each sample is re-embedded on its own, as it would be analysed without the
# atlas, and the cells the atlas assigns to CD8.Q are then marked.
sample_umap_panel <- function(igt_file, sample_code) {
    igt <- readRDS(sprintf("%s/%s", data_path, igt_file))

    # Via @meta.data: see script/Figure1.R for why.
    for (col in c("annotation_level2", "annotation_level1",
                  "cellID", "sample_code")) {
        igt@meta.data[[col]] <- so_orig@meta.data[colnames(igt), col]
    }
    igt <- igt[, !is.na(igt$annotation_level2)]

    igt_sample <- igt[, igt$sample_code == sample_code] %>%
        NormalizeData(verbose = FALSE) %>%
        FindVariableFeatures(verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = 20, verbose = FALSE) %>%
        FindNeighbors(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
        RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE)

    cells_to_highlight <- igt_sample@meta.data %>%
        filter(annotation_level2 == "CD8.Q") %>%
        pull(cellID)

    DimPlot(igt_sample, reduction = "umap",
            cells.highlight = cells_to_highlight,
            cols.highlight = "green", sizes.highlight = 4, pt.size = 1) +
        NoLegend() + NoGrid()
}

panel_pdf(figure_dir, "4f_IGT40_UMAP_CD8Q", 5, 5)
print(sample_umap_panel("IGT40_seurat.Rds", SAMPLE_LCMV))
dev.off()

panel_pdf(figure_dir, "4i_IGT21_UMAP_CD8Q", 5, 5)
print(sample_umap_panel("IGT21_seurat.Rds", SAMPLE_NAIVE))
dev.off()

# ============================================================
# 4g/4j: CD8.Q cells in CD103 vs ITGB7 protein space
# ============================================================
# The whole experiment's cells, with the sample's CD8.Q cells picked out.
protein_panel <- function(igt_id, sample_code) {
    igt <- so_orig[, so_orig$IGT == igt_id]
    igt <- NormalizeData(igt, assay = "ADT",
                         normalization.method = "LogNormalize", verbose = FALSE)
    igt$is_CD8Q_in_sample <-
        igt$annotation_level2 == "CD8.Q" & igt$sample_code == sample_code

    MyFeatureScatter(
        igt, assay = "ADT", slot = "data",
        feature1 = "ITB7", feature2 = "CD103",
        group.by = "is_CD8Q_in_sample",
        highlight_size = 2, highlight_alpha = 1, raster = FALSE,
        cols = ZemmourLib::immgent_colors$level2["CD8.Q"]
    ) + ggtitle(sample_code)
}

panel_pdf(figure_dir, "4g_IGT40_CD8Q_CD103_ITGB7", 5, 5)
print(protein_panel("IGT40", SAMPLE_LCMV))
dev.off()

panel_pdf(figure_dir, "4j_IGT21_CD8Q_CD103_ITGB7", 5, 5)
print(protein_panel("IGT21", SAMPLE_NAIVE))
dev.off()
