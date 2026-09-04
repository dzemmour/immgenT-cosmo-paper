# Extended Data Figure 7. Ablation tests assess the ability of T-RBI to detect
# T-cell states absent from the reference.
#
# Individual CD4 clusters were removed from the reference in turn and their
# cells mapped back as queries, then compared with mapping against the complete
# reference.
#
# Panels produced:
#   S7a  Query cells from each omitted cluster on the ablated CD4 MDE.
#   S7b  Cluster assignments of cells from each omitted cluster.
#   S7c  Annotation confidence for misclassified cells: no ablation, and
#        ablation before and after fine-tuning.
#
# Panel S7d (single-cell discovery scores per omitted cluster) is not
# reproduced here: the cached ablation objects were dieted down to the
# confidence and assignment columns and do not retain the discovery score.
# See analysis/FigureS7.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "TRBI ablation analysis (ok)"
# (line 2642). Three fixes to the Rmd's code were needed:
#   - it read CD4Ablation_trbi_seurat_objects.Rds from out_dir, not data_dir;
#   - it re-read a differently-named atlas object with a bare relative path,
#     silently replacing so_orig for everything below it;
#   - its panel-a loop referenced a `bg_cl` layer built from that object.
# The background here is the CD4 raster from ZemmourLib plus the omitted
# cluster's own position in the complete reference, which is what `bg_cl`
# was for.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds        [primary input]
#   CD4Ablation_trbi_seurat_objects.Rds      [primary input]
#   CD4NoAblation_trbi_seurat.Rds            [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(ggplot2)
    library(scattermore)
    library(tidyr)
    library(scales)
})

source("code/R/setup.R")

figure_dir <- "Extended Data Figure 7"

# The 18 CD4 clusters ablated in turn.
ABLATED_CLUSTERS <- c(
    "CD4.A", "CD4.B", "CD4.C", "CD4.D", "CD4.E", "CD4.F", "CD4.G", "CD4.H",
    "CD4.I", "CD4.Q", "CD4.R", "CD4.S", "CD4.T", "CD4.U", "CD4.X", "CD4.Y",
    "CD4.wZ", "CD4.P"
)

# ============================================================
# Load data
# ============================================================
so_list <- readRDS(sprintf("%s/CD4Ablation_trbi_seurat_objects.Rds",
                           data_path))
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so_cd4 <- so_orig[, so_orig@meta.data %>%
                      filter(annotation_level1 == "CD4",
                             organ_simplified != "thymus") %>%
                      rownames()]

# ============================================================
# S7a: Omitted-cluster cells on the ablated CD4 MDE
# ============================================================
# Two background layers: the whole CD4 reference in black, and in white the
# position the omitted cluster occupied in the complete reference. The query
# cells land in that white gap, which is the panel's point -- MDE preserves
# local neighbourhoods, so cells return to the hole their cluster left behind.
bg_cd4 <- AddEmbeddingRasterBackground(
    lineage = "CD4", interpolate = FALSE, color = "black", alpha = 0.5
)

for (mycl in ABLATED_CLUSTERS) {
    so <- so_list[[mycl]]

    bg_cluster <- geom_scattermore(
        data = as.data.frame(
            so_cd4[, so_cd4$annotation_level2 == mycl][["mde_incremental"]]@cell.embeddings
        ),
        aes(mdeincremental_1, mdeincremental_2),
        color = "white"
    )

    p <- DimPlot(
        so[, so$level2_orig == mycl],
        reduction = "mde_incremental_level2",
        group.by = "level2_orig",
        pt.size = 2, label = FALSE, raster = TRUE,
        cols = mypal_level2
    )
    p$layers <- c(list(bg_cd4), list(bg_cluster), p$layers)

    panel_pdf(figure_dir, sprintf("S7a_ablation_MDE_%s", mycl), 6, 6)
    print(
        p + theme_void() +
            coord_fixed(xlim = c(-2.5, 3.5), ylim = c(-2.5, 2.5)) +
            ggtitle(mycl) + NoLegend()
    )
    dev.off()
}

# ============================================================
# S7b: Cluster assignments after ablation
# ============================================================
# Rows are the omitted cluster, columns the assignment T-RBI gave its cells
# once that cluster was gone. Reported before fine-tuning, which is the
# conservative reading: fine-tuning uses confidently annotated query cells to
# resolve the remainder and so increases reassignment to neighbours.
prop_df <- imap_dfr(so_list, function(so, mycl) {
    so@meta.data %>%
        filter(level2_orig == mycl) %>%
        count(l2 = level2_before_iter, name = "n") %>%
        mutate(level2_orig = mycl, n_total = sum(n), prop = n / n_total)
})

# Columns start with "not classified" (the outcome the ablation is testing for)
# and end with the miniverse clusters, which are not lineage neighbours of
# anything and would otherwise interrupt the block of CD4 clusters.
row_order <- intersect(names(so_list), unique(prop_df$l2))
col_order <- c(row_order, setdiff(unique(prop_df$l2), row_order))
col_order <- c(
    "not classified",
    setdiff(col_order, c("not classified", "CD4.wM", "CD4.wN", "CD4.wZ")),
    c("CD4.wZ", "CD4.wM", "CD4.wN")
)

prop_df <- prop_df %>%
    mutate(
        level2_orig = factor(level2_orig, levels = rev(row_order)),
        l2 = factor(l2, levels = col_order)
    )

p_s7b <- ggplot(prop_df, aes(l2, level2_orig)) +
    geom_point(aes(size = prop, color = prop)) +
    geom_text(aes(label = percent(prop, accuracy = 1)), size = 2.5,
              color = "black") +
    scale_size_continuous(range = c(2, 9), labels = percent,
                          name = "Proportion") +
    scale_color_gradient(low = "white", high = "#B40426", labels = percent,
                         name = "Proportion") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Assignment after ablation", y = "Omitted cluster") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3)) +
    NoGrid()

panel_pdf(figure_dir, "S7b_ablation_confusion_matrix", 8, 6)
print(p_s7b)
dev.off()

# ============================================================
# S7c: Confidence scores for misclassified cells
# ============================================================
# The same cells under three conditions. If two neighbouring clusters were an
# arbitrary split of one state, removing one should reassign its cells to the
# other at undiminished confidence; the drop shown here is the evidence that
# they are distinct.
so_noabl <- readRDS(sprintf("%s/CD4NoAblation_trbi_seurat.Rds", data_path))

df <- lapply(ABLATED_CLUSTERS, function(mycl) {
    so <- so_list[[mycl]]
    cells <- colnames(so)[so$level2_orig == mycl]
    cells <- intersect(cells, colnames(so_noabl))
    md <- so@meta.data[cells, , drop = FALSE]

    # Cells left unclassified have no confidence score to report, so they enter
    # as NA rather than as a low score.
    data.frame(
        cluster = mycl,
        cell = cells,
        no_ablation = so_noabl@meta.data[cells, "confidence_score_after_iter"],
        before_iter = ifelse(
            !is.na(md$level2_before_iter) &
                md$level2_before_iter != "not classified",
            md$confidence_score_before_iter, NA_real_
        ),
        after_iter = ifelse(
            !is.na(md$level2_after_iter) &
                md$level2_after_iter != "not classified",
            md$confidence_score_after_iter, NA_real_
        )
    )
}) %>%
    bind_rows() %>%
    tidyr::pivot_longer(cols = c(no_ablation, before_iter, after_iter),
                        names_to = "condition",
                        values_to = "confidence_score")

df$cluster <- factor(df$cluster, levels = ABLATED_CLUSTERS)
df$condition <- factor(
    df$condition,
    levels = c("no_ablation", "before_iter", "after_iter"),
    labels = c("No ablation", "Ablation - before fine tuning",
               "Ablation - after fine tuning")
)

p_s7c <- ggplot(df, aes(cluster, confidence_score, fill = condition)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
                 outlier.shape = NA) +
    scale_fill_manual(values = c(
        "No ablation" = "white",
        "Ablation - before fine tuning" = "grey70",
        "Ablation - after fine tuning" = "grey30"
    )) +
    labs(x = NULL, y = "Confidence score of misclassified cells", fill = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

panel_pdf(figure_dir, "S7c_ablation_confidence_scores", 10, 6)
print(p_s7c)
dev.off()
