# Figure 8. Transcription factor expression across T-cell states.
#
# Panels produced:
#   8a  Cumulative TF expression across ranked lineages, three most
#       lineage-specific factors labelled.
#   8b  The same across clusters.
#   8c-8m  Expression of Foxp3, Zbtb16, Sox13, Rorc, Pparg, Hes1, Gata3, Tcf7,
#       Tbx21, Tox and Tox2 across clusters.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Calculate pseudobulk (on Midway)"
# (line 2274, which despite its title reads the cached pseudobulk objects) and
# "Select TF expression skyline plots". The Rmd's own header for this material,
# at line 2271, calls it Figure 7; the published order has transcription
# factors as Figure 8.
#
# The Rmd's skyline chunk plotted fifteen factors in one pass: the eleven
# published here plus Sox4, Zbtb32, Tshz2 and Hlf, which are Extended Data
# Figure 10c-f and are drawn by script/FigureS10.R.
#
# As for Figure 7, the Rmd's single Gene_AUC() call per resolution has been
# decomposed into code/R/gene_auc.R; its row-normalised heatmap variant and its
# per-top-gene barplots are not published.
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
})

source("code/R/setup.R")
source("code/R/gene_auc.R")

figure_dir <- "Figure 8"
out_dir <- "output/Figure8"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# ============================================================
# 8a: TF specificity across lineages
# ============================================================
res_l1 <- gene_specificity_auc(pb_level1, tf_all, expr_threshold = 0.5)
write.table(res_l1$auc, file.path(out_dir, "8a_TF_AUC_lineage.txt"),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Only the three most specific factors are labelled: with eight lineages the
# curves are short and more labels overlap them.
panel_pdf(figure_dir, "8a_TF_cumulative_curves_lineage", 7, 7)
print(plot_gene_auc_curves(res_l1, n_top = 3, mypal = mypal))
dev.off()

# ============================================================
# 8b: TF specificity across clusters
# ============================================================
pb_level2$annotation_level1 <- sub("\\..*", "", pb_level2$annotation_level2)

res_l2 <- gene_specificity_auc(pb_level2, tf_all, expr_threshold = 0.5)
write.table(res_l2$auc, file.path(out_dir, "8b_TF_AUC_cluster.txt"),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

panel_pdf(figure_dir, "8b_TF_cumulative_curves_cluster", 7, 7)
print(plot_gene_auc_curves(res_l2, n_top = 20, mypal = mypal))
dev.off()

# ============================================================
# 8c-8m: Selected transcription factors across clusters
# ============================================================
# In panel order.
tf_panels <- c("8c" = "Foxp3", "8d" = "Zbtb16", "8e" = "Sox13",
               "8f" = "Rorc",  "8g" = "Pparg",  "8h" = "Hes1",
               "8i" = "Gata3", "8j" = "Tcf7",   "8k" = "Tbx21",
               "8l" = "Tox",   "8m" = "Tox2")

lineage_of_cluster <- setNames(
    sub("\\..*", "", colnames(pb_level2)),
    colnames(pb_level2)
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
