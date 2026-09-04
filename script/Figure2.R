# Figure 2. immgenT harmonizes T-cell representation across experiments and
# delineates major T-cell lineages.
#
# Panels produced:
#   2a  All-T MDE coloured by T-cell lineage.
#   2b  Surface marker expression by lineage: TCRb vs TCRgd, CD4 vs CD8b,
#       CD8a vs CD8b.
#   2c  All-T MDE coloured by canonical lineage transcripts.
#   2e  Mki67 on the all-T MDE, marking proliferating cells.
#   2f  CD62L/CD44 protein gate and the two single-positive populations on the
#       all-T MDE.
#   2h  Colonic T cells from two independent experiments, before and after
#       integration.
#
# Panels 2d (clonotype sharing between lineages) and 2g (expanded clonotypes)
# are produced by the TCR analysis and are not in scope here.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Figure 2/S3 - Integration (ok)"
# (line 1009) and "Colon in IGT20, IGT24 and IGT27 (Fig 2) (ok)" (line 1113).
# That second header names IGT24, but the code under it only ever handled IGT20
# and IGT27, which are the two experiments the published panel shows.
#
# The same Rmd sections also produce Extended Data Figure 3 panels; those live
# in script/FigureS3.R.
#
# Panel 2c: the Rmd passed the six transcripts in the order Foxp3, Cd4, Cd8b1,
# Cd8a, Trdc, Zbtb16. Reordered here to the published layout (Trdc, Cd8b1, Cd4,
# Zbtb16, Cd8a, Foxp3); same six genes, same settings.
#
# Panel 2e: the Rmd drew Mki67 as one facet of a five-gene FeaturePlot that also
# included Foxp3, Zbtb16, Cdkn2a and Cdkn2b at max.cutoff = "q90". Only Mki67 is
# published, so it is drawn on its own here, at the same cutoff.
#
# Panel 2f: the Rmd additionally gated CD62L-single-positive, CD44-single-
# positive, double-positive and double-negative populations and drew all four on
# the all-T MDE. Only the two single-positive panels are published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   IGT20_seurat.Rds, IGT27_seurat.Rds  [primary inputs]

suppressPackageStartupMessages({
    library(dplyr)
    library(viridis)
})

source("code/R/setup.R")

figure_dir <- "Figure 2"

# Panel 2f subsamples cells for the protein scatter; fixed so the panel is
# reproducible.
set.seed(1)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus") %>% rownames()]

# ============================================================
# 2a: All-T MDE coloured by lineage
# ============================================================
p_2a <- DimPlot(
    object = so, reduction = "mde2_totalvi_20241006",
    group.by = "annotation_level1", raster = TRUE,
    raster.dpi = c(1024, 1024)
) +
    scale_color_manual(values = mypal_level1) +
    NoGrid()

panel_pdf(figure_dir, "2a_MDE_level1", 5, 5)
print(p_2a + NoLegend())
dev.off()

# ============================================================
# 2b: Surface marker expression by lineage
# ============================================================
# Dashed lines mark the positivity thresholds used to reconcile the
# protein-based lineage calls with the transcriptome-based annotation. They are
# drawn for orientation; the annotation itself is not a product of these gates.
so <- NormalizeData(so, assay = "ADT",
                    normalization.method = "LogNormalize", verbose = FALSE)

features_list <- list(
    c("TCRB", "TCRGD"),
    c("CD4", "CD8B"),
    c("CD8A", "CD8B")
)
thresholds <- c("TCRB" = 3, "TCRGD" = 4, "CD4" = 4, "CD8B" = 5, "CD8A" = 4.5)

for (i in seq_along(features_list)) {
    feature1 <- features_list[[i]][1]
    feature2 <- features_list[[i]][2]
    p <- FeatureScatter(
        so, slot = "data", feature1 = feature1, feature2 = feature2,
        group.by = "annotation_level1", raster = TRUE, pt.size = 1,
        plot.cor = FALSE
    ) +
        scale_color_manual(values = mypal_level1) +
        geom_vline(xintercept = thresholds[[feature1]],
                   linetype = "dashed", color = "brown") +
        geom_hline(yintercept = thresholds[[feature2]],
                   linetype = "dashed", color = "brown") +
        NoLegend()

    panel_pdf(figure_dir, sprintf("2b_FeatureScatter_%s_%s", feature1, feature2),
              3, 3)
    print(p)
    dev.off()
}

# ============================================================
# 2c: Canonical lineage transcripts on the all-T MDE
# ============================================================
# RNA in this object is already log-normalised (see load_immgent()); no further
# normalisation is applied. Expression is capped at the 95th percentile so that
# a handful of very high cells do not flatten the rest of the scale.
panel_pdf(figure_dir, "2c_FeaturePlot_lineage_markers", 12, 15)
print(
    FeaturePlot(
        so, reduction = "mde2_totalvi_20241006",
        features = c("Trdc", "Cd8b1", "Cd4", "Zbtb16", "Cd8a", "Foxp3"),
        order = TRUE, cols = c("lightgrey", "blue"), raster = TRUE,
        max.cutoff = "q95"
    ) +
        ggplot2::xlim(-2.5, 2.5) + ggplot2::ylim(-2.5, 2.5)
)
dev.off()

# ============================================================
# 2e: Mki67 marks proliferating cells
# ============================================================
panel_pdf(figure_dir, "2e_FeaturePlot_Mki67", 5, 5)
print(
    FeaturePlot(
        so, reduction = "mde2_totalvi_20241006", features = "Mki67",
        order = TRUE, cols = c("lightgrey", "blue"), raster = TRUE,
        max.cutoff = "q90"
    ) +
        ggplot2::xlim(-2.5, 2.5) + ggplot2::ylim(-2.5, 2.5)
)
dev.off()

# ============================================================
# 2f: CD62L/CD44 gate and the single-positive populations
# ============================================================
# Restricted to CITE-seq cells, and re-derived from so_orig rather than reusing
# the object above, because the ADT normalisation must be computed on exactly
# the cells being gated for the thresholds to mean the same thing.
so_cite <- so_orig[, so_orig@meta.data %>%
                       filter(organ_simplified != "thymus", cite_seq == TRUE) %>%
                       rownames()]
so_cite <- NormalizeData(so_cite, assay = "ADT",
                         normalization.method = "LogNormalize", verbose = FALSE)

# Thresholds read off the bimodal protein distributions shown in the gating
# panel below.
cd62l_gate <- 4.5
cd44_gate <- 5

adt <- so_cite[["ADT"]]$data
so_cite$is_CD62L_sp <- adt["CD62L", ] > cd62l_gate & adt["CD44", ] < cd44_gate
so_cite$is_CD44_sp  <- adt["CD62L", ] < cd62l_gate & adt["CD44", ] > cd44_gate

# The gating plot itself is drawn on a 50,000-cell subsample: at 500,000+ cells
# the density kernel saturates and the two populations stop being separable by
# eye.
so_cite$sub <- FALSE
so_cite$sub[sample(colnames(so_cite), size = 50000, replace = FALSE)] <- TRUE

panel_pdf(figure_dir, "2f_gate_CD62L_CD44", 5, 5)
print(
    MyFeatureScatter(
        so = so_cite[, so_cite$sub], assay = "ADT", slot = "data",
        feature1 = "CD62L", feature2 = "CD44", group.by = "sub",
        raster = TRUE, highlight_size = 1
    ) +
        geom_vline(xintercept = cd62l_gate, linetype = "dashed",
                   color = "brown") +
        geom_hline(yintercept = cd44_gate, linetype = "dashed",
                   color = "brown") +
        NoLegend()
)
dev.off()

panel_pdf(figure_dir, "2f_MDE_CD62Lpos_CD44neg", 5, 5)
print(
    MyDimPlotHighlightDensity(
        seurat_object = so_cite, umap_to_plot = "mde2_totalvi_20241006",
        group.by = "is_CD62L_sp", raster = TRUE, cols = plasma(10),
        highlight_size = 0.1, highlight_pixels = c(512, 512)
    ) + ggtitle("CD44- CD62L+")
)
dev.off()

panel_pdf(figure_dir, "2f_MDE_CD44pos_CD62Lneg", 5, 5)
print(
    MyDimPlotHighlightDensity(
        seurat_object = so_cite, umap_to_plot = "mde2_totalvi_20241006",
        group.by = "is_CD44_sp", raster = TRUE, cols = plasma(10),
        highlight_size = 0.1, highlight_pixels = c(512, 512)
    ) + ggtitle("CD44+ CD62L-")
)
dev.off()

rm(so_cite)

# ============================================================
# 2h: Colonic T cells before and after integration
# ============================================================
# The point of the panel is that two experiments which cluster colonic T cells
# differently on their own land on the same states once integrated. So the
# left-hand plot uses each experiment's own UMAP and its own de novo cluster
# numbering, and the right-hand plot places the same cells on the shared all-T
# MDE coloured by immgenT lineage.
colon_experiments <- list(
    IGT20 = "Colon",       # sample_name pattern selecting the colonic samples
    IGT27 = "No_colon_M"
)

for (igt_id in names(colon_experiments)) {
    sample_pattern <- colon_experiments[[igt_id]]
    igt <- readRDS(sprintf("%s/%s_seurat.Rds", data_path, igt_id))

    # See script/Figure1.R for why these assignments go through @meta.data.
    igt@meta.data$annotation_level1 <- as.character(
        so_orig$annotation_level1[match(colnames(igt), colnames(so_orig))]
    )
    colon_cells <- colnames(igt)[grepl(sample_pattern, igt$sample_name)]

    # A palette keyed on this experiment's own cluster labels, so the same
    # cluster keeps its colour between the pre- and post-integration plots.
    mypal_batch <- setNames(mypal, unique(igt$RNA_clusters))
    mypal_batch <- mypal_batch[!is.na(names(mypal_batch))]

    panel_pdf(figure_dir,
              sprintf("2h_%s_UMAP_RNAclusters", igt_id), 5, 5)
    print(
        DimPlot(igt[, colon_cells], reduction = "umap_rna",
                group.by = "RNA_clusters") +
            scale_color_manual(values = mypal_batch) + NoLegend() + NoGrid()
    )
    dev.off()

    # Carry this experiment's cluster labels onto the atlas object so the same
    # cells can be highlighted on the shared MDE.
    tmp <- so_orig
    tmp@meta.data$RNA_clusters <- as.character(
        igt$RNA_clusters[match(colnames(tmp), colnames(igt))]
    )

    panel_pdf(figure_dir,
              sprintf("2h_%s_MDE_RNAclusters", igt_id), 5, 5)
    print(
        highlight_mde(
            tmp, umap_to_plot = "mde2_totalvi_20241006",
            cells_to_highlight = colon_cells,
            highlight_column_name = "RNA_clusters", mycols = mypal_batch,
            labelclusters = FALSE, which = "plot2"
        )
    )
    dev.off()

    panel_pdf(figure_dir,
              sprintf("2h_%s_MDE_level1", igt_id), 5, 5)
    print(
        highlight_mde(
            tmp, umap_to_plot = "mde2_totalvi_20241006",
            cells_to_highlight = colon_cells,
            highlight_column_name = "annotation_level1", mycols = mypal_level1,
            labelclusters = FALSE, which = "plot2"
        )
    )
    dev.off()

    rm(tmp, igt)
}
