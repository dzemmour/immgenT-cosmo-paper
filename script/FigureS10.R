# Extended Data Figure 10. Transcription factor expression in immgenT clusters.
#
# Panels produced:
#   S10a  Heatmap of the 20 most lineage-specific transcription factors.
#   S10b  Heatmap of the 20 most cluster-specific transcription factors.
#   S10c  Sox4 across clusters.
#   S10d  Zbtb32 across clusters.
#   S10e  Tshz2 across clusters.
#   S10f  Hlf across clusters.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Calculate pseudobulk (on Midway)"
# (line 2274) and its "Select TF expression skyline plots" chunk, which drove
# both this figure and Figure 8. The heatmaps here are two of the three variants
# the Rmd's Gene_AUC() emitted per resolution; the row-normalised variant is not
# published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                            [primary input]
#   igt1_96_pseudobulk_byclusterannotationlevel1_log1pCP10K.Rds  [primary input]
#   igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds  [primary input]
#   TF_list.txt                                                  [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
})

source("code/R/setup.R")
source("code/R/gene_auc.R")

figure_dir <- "Extended Data Figure 10"

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
tf_all <- read.table(sprintf("%s/TF_list.txt", data_path),
                     header = FALSE, stringsAsFactors = FALSE)[, 1]

pb_level1 <- readRDS(sprintf(
    "%s/igt1_96_pseudobulk_byclusterannotationlevel1_log1pCP10K.Rds", data_path
))
pb_level2 <- readRDS(sprintf(
    "%s/igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds", data_path
))

for (pb_name in c("pb_level1", "pb_level2")) {
    pb <- get(pb_name)
    if (!"annotation_level1" %in% colnames(pb@meta.data)) {
        pb$annotation_level1 <- sub("\\..*", "", pb$annotation_level2)
    }
    pb$annotation_level2_group <- as.character(
        so_orig$annotation_level2_group[
            match(colnames(pb), so_orig$annotation_level2)
        ]
    )
    assign(pb_name, pb)
}

ann_palettes <- list(
    annotation_level1 = mypal_level1,
    annotation_level2_group = mypal_level2group
)

# ============================================================
# S10a: Most lineage-specific transcription factors
# ============================================================
# Columns fixed in lineage order rather than clustered: with eight columns the
# useful comparison is between named lineages, not between dendrogram
# neighbours.
res_l1 <- gene_specificity_auc(pb_level1, tf_all, expr_threshold = 0.5)

ht_s10a <- plot_gene_heatmap(
    res_l1, pb_level1, n_top = 20,
    col_order = levels(so_orig$annotation_level1),
    colours = ColorRamp,
    ann_cols = "annotation_level1",
    ann_palettes = ann_palettes,
    main = "20 most lineage-specific TFs, log1p CP10K"
)

panel_pdf(figure_dir, "S10a_heatmap_TF_top20_lineage", 7, 7)
grid::grid.newpage()
grid::grid.draw(ht_s10a$gtable)
dev.off()

# ============================================================
# S10b: Most cluster-specific transcription factors
# ============================================================
res_l2 <- gene_specificity_auc(pb_level2, tf_all, expr_threshold = 0.5)

ht_s10b <- plot_gene_heatmap(
    res_l2, pb_level2, n_top = 20,
    col_order = levels(so_orig$annotation_level2),
    colours = colorRampPalette(c("white", "orange", "brown"))(100),
    ann_cols = c("annotation_level1", "annotation_level2_group"),
    ann_palettes = ann_palettes,
    main = "20 most cluster-specific TFs, log1p CP10K"
)

panel_pdf(figure_dir, "S10b_heatmap_TF_top20_cluster", 12, 7)
grid::grid.newpage()
grid::grid.draw(ht_s10b$gtable)
dev.off()

# ============================================================
# S10c-S10f: Selected transcription factors across clusters
# ============================================================
tf_panels <- c("S10c" = "Sox4", "S10d" = "Zbtb32",
               "S10e" = "Tshz2", "S10f" = "Hlf")

lineage_of_cluster <- setNames(
    sub("\\..*", "", colnames(pb_level2)), colnames(pb_level2)
)
cluster_levels <- intersect(levels(so_orig$annotation_level2),
                            colnames(pb_level2))

for (panel in names(tf_panels)) {
    gene <- tf_panels[[panel]]
    panel_pdf(figure_dir, sprintf("%s_barplot_%s", panel, gene), 12, 5)
    print(plot_gene_barplot(
        pb_level2, gene,
        group_levels = cluster_levels,
        group_colors = mypal_level1,
        fill_by = lineage_of_cluster
    ))
    dev.off()
}
