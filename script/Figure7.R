# Figure 7. Effector gene expression across T-cell states.
#
# Panels produced:
#   7a  Heatmap of curated effector gene expression across clusters, rows and
#       columns hierarchically clustered, annotated by lineage and activation
#       state.
#   7b  Il4 across clusters.
#   7c  Il17a across clusters.
#   7d  Ifng across clusters.
#
# Panels 7e-7g (sample-to-sample variability of Il4 in CD4.T, Il17a in CD4.X
# and Ifng in CD4.Q) require a per-sample pseudobulk matrix that is not among
# the cached inputs; see analysis/Figure7.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 8: Effector molecules"
# (line 2370) and its "Across Level2" subsection. The Rmd's header numbering is
# from an earlier draft in which effector molecules were Figure 8 and
# transcription factors Figure 7; the published order is the reverse.
#
# The Rmd drove both this figure and Extended Data Figure 9a through one
# Gene_AUC() call per pseudobulk resolution, which printed the cumulative
# curves, a barplot per top gene and three heatmap variants in a single pass.
# That function has been decomposed into code/R/gene_auc.R so each published
# panel is produced by its own explicit call; the effector-gene cumulative
# curves and the row-normalised heatmap variant it also emitted are not
# published.
#
# The bar plots 7b-7d were not in the Rmd as standalone calls -- Gene_AUC()
# emitted a barplot for each of its top-scoring genes, which for the effector
# list does not put Il4, Il17a and Ifng in the top 20. They are drawn here
# directly with the same helper Figure 8 uses for its named transcription
# factors, which is the same plotting code Gene_AUC() used internally.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                            [primary input]
#   igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds  [primary input]
#   effmol_curated.txt                                           [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
})

source("code/R/setup.R")
source("code/R/gene_auc.R")

figure_dir <- "Figure 7"
out_dir <- "output/Figure7"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
effmol <- read.table(sprintf("%s/effmol_curated.txt", data_path),
                     header = TRUE)[, 1]

pseudobulk <- readRDS(sprintf(
    "%s/igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds", data_path
))

# The pseudobulk object carries only the cluster name, so lineage and
# activation state are recovered here: lineage from the cluster name (every
# cluster is "<lineage>.<letter>"), activation state by lookup against the
# atlas.
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

# Fill colour per cluster, taken from that cluster's lineage, so the bar plots
# below read as lineage blocks.
lineage_of_cluster <- setNames(
    sub("\\..*", "", colnames(pseudobulk)),
    colnames(pseudobulk)
)

# ============================================================
# 7a: Effector gene expression heatmap
# ============================================================
# The full curated effector list, less the genes that never reach 0.5 log1p
# CP10K in any cluster. Hierarchically clustered on both axes: the point of the
# panel is that effector genes group into modules reused across lineages, which
# a fixed cluster order would hide.
#
# Extended Data Figure 9a is the same plot over the 1,821-gene GO:0005615
# (extracellular space) list, cut to its 100 most cluster-specific genes.
res <- gene_specificity_auc(pseudobulk, effmol, expr_threshold = 0.5)
message("effector genes passing the expression filter: ", nrow(res$auc))
write.table(res$auc, file.path(out_dir, "7a_effector_gene_AUC.txt"),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

ht_7a <- plot_gene_heatmap(
    res, pseudobulk, n_top = nrow(res$auc),
    colours = colorRampPalette(c("white", "red", "brown"))(100),
    ann_cols = c("annotation_level1", "annotation_level2_group"),
    ann_palettes = ann_palettes,
    main = "Effector genes, expression log1p CP10K"
)

panel_pdf(figure_dir, "7a_heatmap_effector_genes", 15, 15)
grid::grid.newpage()
grid::grid.draw(ht_7a$gtable)
dev.off()

# ============================================================
# 7b-7d: Individual effector cytokines across clusters
# ============================================================
# Three cytokines whose textbook assignment is to a single T-helper subset.
# Drawn across every cluster to show that each is instead produced by clusters
# from several lineages.
cytokines <- c("7b" = "Il4", "7c" = "Il17a", "7d" = "Ifng")

for (panel in names(cytokines)) {
    gene <- cytokines[[panel]]
    panel_pdf(figure_dir, sprintf("%s_barplot_%s", panel, gene), 12, 5)
    print(plot_gene_barplot(
        pseudobulk, gene,
        group_levels = intersect(levels(so_orig$annotation_level2),
                                 colnames(pseudobulk)),
        group_colors = mypal_level1,
        fill_by = lineage_of_cluster
    ))
    dev.off()
}
