# Extended Data Figure 9. Effector molecule expression in immgenT clusters.
#
# Panels produced:
#   S9a  Heatmap of the 100 most cluster-specific genes from the GO:0005615
#        (extracellular space) list, across clusters.
#   S9b  Il9 across clusters.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 8: Effector molecules"
# (line 2370), which drove this figure and Figure 7 through the same Gene_AUC()
# calls. Figure 7a uses the 59-gene curated effector list; this panel uses the
# 1,821-gene GO:0005615 list, cut to its 100 most cluster-specific genes, which
# is what distinguishes the two heatmaps.
#
# Il9 is in the curated list but not among its top-scoring genes, so the Rmd's
# Gene_AUC() barplot sweep never emitted it; it is drawn directly here.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                            [primary input]
#   igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds  [primary input]
#   effmol_GO0005615_curated.txt                                 [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
})

source("code/R/setup.R")
source("code/R/gene_auc.R")

figure_dir <- "Extended Data Figure 9"
out_dir <- "output/FigureS9"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
effmol_go <- read.table(sprintf("%s/effmol_GO0005615_curated.txt", data_path),
                        header = TRUE)[, 1]

pseudobulk <- readRDS(sprintf(
    "%s/igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds", data_path
))
pseudobulk$annotation_level1 <- sub("\\..*", "", pseudobulk$annotation_level2)
pseudobulk$annotation_level2_group <- as.character(
    so_orig$annotation_level2_group[
        match(pseudobulk$annotation_level2, so_orig$annotation_level2)
    ]
)

ann_palettes <- list(
    annotation_level1 = mypal_level1,
    annotation_level2_group = mypal_level2group
)

# ============================================================
# S9a: Top 100 cluster-specific extracellular genes
# ============================================================
res <- gene_specificity_auc(pseudobulk, effmol_go, expr_threshold = 0.5)
message("GO:0005615 genes passing the expression filter: ", nrow(res$auc))
write.table(res$auc, file.path(out_dir, "S9a_effector_GO0005615_AUC.txt"),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

ht_s9a <- plot_gene_heatmap(
    res, pseudobulk, n_top = 100,
    colours = colorRampPalette(c("white", "red", "brown"))(100),
    ann_cols = c("annotation_level1", "annotation_level2_group"),
    ann_palettes = ann_palettes,
    main = "GO:0005615 genes, expression log1p CP10K"
)

panel_pdf(figure_dir, "S9a_heatmap_effector_GO0005615_top100", 15, 15)
grid::grid.newpage()
grid::grid.draw(ht_s9a$gtable)
dev.off()

# ============================================================
# S9b: Il9 across clusters
# ============================================================
lineage_of_cluster <- setNames(
    sub("\\..*", "", colnames(pseudobulk)), colnames(pseudobulk)
)

panel_pdf(figure_dir, "S9b_barplot_Il9", 12, 5)
print(plot_gene_barplot(
    pseudobulk, "Il9",
    group_levels = intersect(levels(so_orig$annotation_level2),
                             colnames(pseudobulk)),
    group_colors = mypal_level1,
    fill_by = lineage_of_cluster
))
dev.off()
