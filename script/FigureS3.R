# Extended Data Figure 3. Spleen standards across experiments; T-cell lineage
# separation.
#
# Panels produced:
#   S3a  All-T MDE highlighting the spleen standard samples, coloured by
#        experiment.
#   S3b  All-T MDE highlighting each lineage separately.
#   S3c  Surface marker expression, one column per lineage.
#   S3d  iNKT cells, identified by their invariant TCR.
#   S3e  MAIT cells, identified by their invariant TCR.
#   S3f  Proliferating and miniverse cells on the all-T MDE.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 2/S3 - Integration (ok)"
# (line 1009), whose chunks produce both Figure 2 and this figure.
#
# Panels S3a and S3b highlight tens to hundreds of thousands of cells, so their
# highlight layers are rasterised (highlight_raster = TRUE, at the same 1024 px
# as the background) as every other highlight panel in this repository does.
# Left as vector points, the CD4 and CD8 panels alone came to 13 MB and 11 MB of
# overlapping paths -- files no vector editor will open, for no gain, since the
# point cloud carries no structure worth editing. Titles, legends and axes stay
# vector.
#
# Panel S3f: the Rmd drew all five annotation_level2_group values in one
# DimPlot with the resting/activated palette. The published panel shows only the
# proliferating and miniverse groups, in black and red, so it is drawn here as a
# two-group highlight over the rest of the atlas in grey.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   color_palette_igt.csv               [curated input]

suppressPackageStartupMessages({
    library(dplyr)
})

source("code/R/setup.R")

figure_dir <- "Extended Data Figure 3"

LINEAGES <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus") %>% rownames()]
mypal_igt <- load_igt_palette()

# ============================================================
# S3a: Spleen standard samples across experiments
# ============================================================
# A standardised spleen sample was spiked into every experiment; those cells
# are highlighted here and coloured by experiment.
panel_pdf(figure_dir, "S3a_MDE_spleen_standard_by_IGT", 5, 5)
print(highlight_mde(
    seurat_object = so, umap_to_plot = "mde2_totalvi_20241006",
    cells_to_highlight = names(which(so$spleen_standard == TRUE)),
    highlight_column_name = "IGT", pixels = c(1024, 1024),
    mycols = mypal_igt, title = "Spleen standards",
    highlight_size = 1, highlight_alpha = 1, labelclusters = FALSE,
    highlight_raster = TRUE, highlight_pixels = c(1024, 1024),
    which = "plot1"
))
dev.off()

# ============================================================
# S3b: Each lineage on the all-T MDE
# ============================================================
# One panel per lineage rather than the single colour-coded MDE of Figure 2a,
# at full opacity, with the lineage name as a title rather than as labels
# inside the panel.
for (lin in LINEAGES) {
    panel_pdf(figure_dir, sprintf("S3b_MDE_%s", lin), 10, 10)
    print(
        highlight_mde(
            seurat_object = so, umap_to_plot = "mde2_totalvi_20241006",
            cells_to_highlight = colnames(so)[so$annotation_level1 == lin],
            highlight_column_name = "annotation_level1",
            pixels = c(1024, 1024),
            mycols = mypal_level1,
            highlight_size = 1, highlight_alpha = 1, labelclusters = FALSE,
            highlight_raster = TRUE, highlight_pixels = c(1024, 1024),
            which = "plot2"
        ) + ggtitle(lin)
    )
    dev.off()
}

# ============================================================
# S3c: Surface marker expression per lineage
# ============================================================
# The same three marker pairs as Figure 2b, split by lineage.
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
        group.by = "annotation_level1", split.by = "annotation_level1",
        raster = TRUE, pt.size = 1
    ) +
        scale_color_manual(values = mypal_level1) +
        geom_vline(xintercept = thresholds[[feature1]],
                   linetype = "dashed", color = "brown") +
        geom_hline(yintercept = thresholds[[feature2]],
                   linetype = "dashed", color = "brown") +
        NoLegend()

    panel_pdf(figure_dir,
              sprintf("S3c_FeatureScatter_%s_%s_bylineage", feature1, feature2),
              15, 3)
    print(p)
    dev.off()
}

# ============================================================
# S3d/S3e: Invariant-TCR populations
# ============================================================
# Both populations are called from TCR sequence alone, with no reference to the
# transcriptome; the selection criteria are in the plot titles.
panel_pdf(figure_dir, "S3d_MDE_iNKT", 10, 10)
print(highlight_mde(
    seurat_object = so, umap_to_plot = "mde2_totalvi_20241006",
    cells_to_highlight = so@meta.data %>% filter(iNKT) %>% rownames(),
    highlight_column_name = "iNKT", mycols = c("grey", "red"),
    labelclusters = FALSE, highlight_raster = TRUE,
    title = paste("iNKT: TRAV11 + TRAJ18 + CDR3a junction",
                  "CVVGDRGSALGRLHF / CVVADRGSALGRLHF / CVVVDRGSALGRLHF"),
    which = "plot1"
))
dev.off()

panel_pdf(figure_dir, "S3e_MDE_MAIT", 10, 10)
print(highlight_mde(
    seurat_object = so, umap_to_plot = "mde2_totalvi_20241006",
    cells_to_highlight = so@meta.data %>% filter(MAIT) %>% rownames(),
    highlight_column_name = "MAIT", mycols = c("grey", "blue"),
    labelclusters = FALSE, highlight_raster = TRUE,
    title = "MAIT: TRAV1 + TRAJ33 + CDR3a length 12 (beta chain indifferent)",
    which = "plot1"
))
dev.off()

# ============================================================
# S3f: Proliferating and miniverse cells
# ============================================================
# The two annotation_level2_group values shown by this panel, highlighted
# together over the rest of the atlas in grey.
so$proliferating_or_miniverse <- dplyr::case_when(
    so$annotation_level2_group == "proliferating" ~ "proliferating",
    so$annotation_level2_group == "miniverse" ~ "miniverse",
    TRUE ~ NA_character_
)

panel_pdf(figure_dir, "S3f_MDE_proliferating_miniverse", 5, 5)
print(highlight_mde(
    seurat_object = so, umap_to_plot = "mde2_totalvi_20241006",
    cells_to_highlight =
        colnames(so)[!is.na(so$proliferating_or_miniverse)],
    highlight_column_name = "proliferating_or_miniverse",
    pixels = c(1024, 1024),
    mycols = c(proliferating = "black", miniverse = "red"),
    title = "Proliferating (black) and miniverse (red)",
    highlight_size = 0.3, highlight_alpha = 1,
    highlight_raster = TRUE, labelclusters = FALSE,
    which = "plot1"
))
dev.off()
