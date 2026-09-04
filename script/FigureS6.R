# Extended Data Figure 6. T-RBI metrics.
#
# Panels produced:
#   S6b  Non-T cells from the external studies on the all-T MDE, as a positive
#        control for detecting states absent from the reference.
#   S6c  Proportion of cells per study without a lineage (level-1) annotation.
#   S6d  Lineage marker expression across the integrated external cells.
#   S6e  Proportion of T cells per study without a cluster (level-2) annotation.
#   S6f  Cluster-level unannotated T cells on the all-T MDE.
#   S6g  Single-cell discovery score on the all-T MDE.
#   S6h  Discovery score distribution per study.
#
# Panel S6a (UMAP of the Miller et al. dataset with unclassified cells marked)
# needs that study's own embedding, which is not among the cached inputs; see
# analysis/FigureS6.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 5 T-RBI" (line 2077), chunks
# "FeaturePlots", "showing nonT outside the MDE", "Barplot with % annotated",
# the not-classified MDE chunk and "Discovery score of not classified cells".
# The Rmd's not-classified chunk also drew the level-1 version of panel f,
# which is not published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                              [primary input]
#   trbi_17studies_diet_merged.Rds                                 [primary input]
#   TRBI_discovery_scores_tbl_merged_DatasetPCA_withnonT.Rds       [code/pipeline/03]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scattermore)
})

source("code/R/setup.R")
source("code/R/trbi_plots.R")

figure_dir <- "Extended Data Figure 6"

# level2_final values that name a lineage rather than a cluster: cells placed
# confidently in a lineage but left unresolved at cluster level.
UNANNOTATED_LEVEL2 <- c("not classified", "CD4", "CD8", "Treg", "gdT",
                        "nonconv", "Tz", "CD8aa", "DN", "DP", "thymocyte")

# ============================================================
# Load data
# ============================================================
so_atlas <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so_merged_orig <- readRDS(sprintf("%s/trbi_17studies_diet_merged.Rds",
                                  data_path))

so_merged <- so_merged_orig[, !so_merged_orig$level1_final %in%
                                c("nonT", "unclear")]
so_merged$level1_final[so_merged$level1_final == "nonconv"] <- "Tz"

bkrg <- trbi_background(so_atlas, "mde2_totalvi_20241006")

# ============================================================
# S6b: Non-T cells as a positive control
# ============================================================
# Contaminating non-T cells are the one population certain to be absent from a
# T-cell reference, so where they land tests whether integration forces query
# cells onto the reference manifold. They do not overlap it.
so_merged_orig$is_nonT <- so_merged_orig$level1_final == "nonT"
dat_all <- trbi_foreground_df(so_merged_orig, "mde_incremental_allT")

panel_pdf(figure_dir, "S6b_MDE_nonT_cells", 5, 5)
print(
    bkrg +
        geom_scattermore(data = dat_all, aes(dim1, dim2), color = "black",
                         pointsize = 3, pixels = c(1024, 1024)) +
        geom_point(data = dat_all[dat_all$is_nonT, ], aes(dim1, dim2),
                   color = "red", size = 1, alpha = 0.1) +
        ggtitle("non-T in red, other query cells black, immgenT grey") +
        trbi_theme() + NoLegend()
)
dev.off()

# ============================================================
# S6c/S6e: Unannotated fractions per study
# ============================================================
prop_level1 <- so_merged@meta.data %>%
    as.data.frame() %>%
    count(dataset, level1_final, name = "n_cells") %>%
    group_by(dataset) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    filter(level1_final == "not classified")

panel_pdf(figure_dir, "S6c_barplot_unannotated_level1", 5, 5)
print(
    ggplot(prop_level1, aes(x = dataset, y = prop * 100)) +
        geom_col(fill = "black") +
        labs(x = "Dataset", y = "% not annotated",
             title = "Cells without a lineage annotation") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        NoGrid()
)
dev.off()

prop_level2 <- so_merged@meta.data %>%
    as.data.frame() %>%
    group_by(dataset) %>%
    summarise(prop_unannotated = mean(level2_final %in% UNANNOTATED_LEVEL2),
              .groups = "drop")

panel_pdf(figure_dir, "S6e_barplot_unannotated_level2", 5, 5)
print(
    ggplot(prop_level2, aes(x = dataset, y = prop_unannotated * 100)) +
        geom_col(fill = "black") +
        labs(x = "Dataset", y = "% not annotated",
             title = "T cells without a cluster annotation") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        NoGrid()
)
dev.off()

# ============================================================
# S6d: Lineage markers across the integrated external cells
# ============================================================
# Marker expression is not used by the mapping, so agreement between where a
# query cell lands and which lineage transcript it expresses is an independent
# check on the annotation.
for (gene in c("Cd4", "Cd8b1", "Cd8a", "Trdc", "Foxp3", "Zbtb16")) {
    # Cd8b1 and Cd8a are left uncapped, as in the original analysis; the other
    # four are capped at the 90th percentile.
    max_cut <- if (gene %in% c("Cd8b1", "Cd8a")) NA else "q90"
    panel_pdf(figure_dir, sprintf("S6d_FeaturePlot_%s", gene), 5, 5)
    print(
        FeaturePlot(so_merged_orig, raster.dpi = c(1024, 1024),
                    features = gene, order = TRUE,
                    reduction = "mde_incremental_allT",
                    max.cutoff = max_cut) +
            xlim(TRBI_XLIM) + ylim(TRBI_YLIM) + theme_void()
    )
    dev.off()
}

# ============================================================
# S6f: Cluster-level unannotated T cells
# ============================================================
# Unannotated T cells are interspersed among annotated ones rather than
# forming a population of their own, which is what distinguishes ordinary
# boundary ambiguity from a state the reference lacks.
so_merged$is_unannotated_level2 <-
    so_merged$level2_final %in% UNANNOTATED_LEVEL2
dat_frg <- trbi_foreground_df(so_merged, "mde_incremental_allT")

panel_pdf(figure_dir, "S6f_MDE_unannotated_level2", 5, 5)
print(
    bkrg +
        geom_scattermore(data = dat_frg, aes(dim1, dim2), col = "grey0",
                         pointsize = 3, pixels = c(1024, 1024)) +
        geom_scattermore(data = dat_frg[dat_frg$is_unannotated_level2, ],
                         aes(dim1, dim2), color = "red",
                         pixels = c(512, 512)) +
        ggtitle("Cluster-level unannotated T cells in red") +
        trbi_theme() + NoLegend()
)
dev.off()

# ============================================================
# S6g/S6h: Discovery score
# ============================================================
# Scores are read from cache: each is computed in its own study's PCA space,
# independently of the integration (see code/R/discovery.R). Non-T cells are
# retained here as the positive control, which is what distinguishes these two
# panels from Figure 5f.
scores <- readRDS(sprintf(
    "%s/TRBI_discovery_scores_tbl_merged_DatasetPCA_withnonT.Rds", data_path
))

# Diverging ramp centred on 1, the value at which a cell is equidistant from
# annotated and unannotated neighbours, and squished at 2 so the conservative
# 1.1 threshold sits in the visible part of the scale.
score_ramp <- rev(colorRampPalette(c("red", "white", "blue"))(20))

panel_pdf(figure_dir, "S6g_MDE_discovery_score", 5, 5)
print(
    bkrg +
        geom_scattermore(data = dat_frg, aes(dim1, dim2), col = "grey0",
                         pointsize = 3, pixels = c(1024, 1024)) +
        geom_point(data = scores[!is.na(scores$ratio), ],
                   aes(mde_incremental_allT_dim1, mde_incremental_allT_dim2,
                       colour = ratio),
                   size = 0.25) +
        scale_color_gradientn(colours = score_ramp, limits = c(0, 2),
                              oob = scales::squish) +
        ggtitle("Discovery score") +
        trbi_theme()
)
dev.off()

panel_pdf(figure_dir, "S6h_discovery_score_by_dataset", 5, 5)
print(
    ggplot(scores, aes(x = dataset, y = ratio)) +
        geom_boxplot(outlier.shape = NA, fill = NA, colour = "black") +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.1) +
        geom_hline(yintercept = 1.1, colour = "brown", linetype = "dashed",
                   linewidth = 0.5) +
        labs(x = "Dataset", y = "Discovery score") +
        ylim(0, 10) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        NoGrid()
)
dev.off()
