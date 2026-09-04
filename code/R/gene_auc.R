# Gene specificity across lineages or clusters, and the plots built from it.
#
# The specificity measure behind Figures 7-8 and Extended Data Figures 9-10.
# For one gene, its pseudobulk expression across groups is scaled to sum to 1
# and the groups are ranked from highest to lowest. The cumulative curve of
# that scaled expression rises steeply for a gene confined to one or two groups
# and follows the diagonal for a gene expressed everywhere; the area under it
# is therefore a specificity score. It is a ROC-like construction, hence the
# name, but computed on ranked expression shares rather than on labels.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
})

#' Specificity score per gene
#'
#' @param pseudobulk Seurat object of pseudobulk profiles, one column per
#'   group, `data` layer holding log1p CP10K values.
#' @param genes Genes to score; those absent from `pseudobulk` are dropped.
#' @param expr_threshold Minimum log1p CP10K a gene must reach in at least one
#'   group to be scored. Removes genes whose apparent specificity comes only
#'   from being near-zero everywhere.
#' @return List with `auc` (one row per gene, descending score), `data_long`
#'   (ranked cumulative curves), `mat` (raw expression) and `mat_norm`
#'   (row-normalised expression).
gene_specificity_auc <- function(pseudobulk, genes, expr_threshold = 0.5) {
    stopifnot("RNA" %in% names(pseudobulk))
    genes <- intersect(genes, rownames(pseudobulk))

    mat <- pseudobulk[["RNA"]]$data[genes, , drop = FALSE]

    # Scale each gene to sum to 1 across groups, so the score reflects how
    # expression is distributed rather than how highly the gene is expressed.
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- NA
    mat_norm <- mat / row_sums

    filt_expr <- rownames(mat)[rowSums(mat > expr_threshold, na.rm = TRUE) >= 1]
    mat_norm_filt <- mat_norm[filt_expr, , drop = FALSE]

    data_long <- mat_norm_filt |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "Gene") |>
        tidyr::pivot_longer(cols = -Gene, names_to = "Group",
                            values_to = "Expression") |>
        dplyr::group_by(Gene) |>
        dplyr::arrange(dplyr::desc(Expression), .by_group = TRUE) |>
        dplyr::mutate(
            CumulativeExpression = cumsum(Expression),
            Rank = dplyr::row_number()
        ) |>
        dplyr::ungroup()

    gene_auc <- data_long |>
        dplyr::group_by(Gene) |>
        dplyr::summarise(auc = sum(CumulativeExpression) / dplyr::n(),
                         .groups = "drop") |>
        dplyr::arrange(dplyr::desc(auc))

    list(auc = gene_auc, data_long = data_long,
         mat = mat, mat_norm = mat_norm)
}

#' Cumulative specificity curves, top genes highlighted
#'
#' All scored genes are drawn as faint grey lines; the `n_top` most specific
#' are drawn in colour and labelled.
#'
#' @param res Output of `gene_specificity_auc()`.
#' @param n_top Number of genes to highlight and label.
#' @param mypal Colours for the highlighted genes; recycled if unnamed.
plot_gene_auc_curves <- function(res, n_top = 20, mypal = NULL) {
    top_genes <- head(res$auc$Gene, n_top)

    cols <- if (is.null(mypal)) {
        setNames(grDevices::rainbow(length(top_genes)), top_genes)
    } else if (is.null(names(mypal))) {
        setNames(rep_len(mypal, length(top_genes)), top_genes)
    } else {
        mypal
    }

    ggplot(res$data_long,
           aes(x = Rank, y = CumulativeExpression, group = Gene)) +
        geom_line(linewidth = 0.5, alpha = 0.1) +
        geom_line(
            data = subset(res$data_long, Gene %in% top_genes),
            aes(colour = Gene), linewidth = 1, alpha = 1
        ) +
        ggrepel::geom_label_repel(
            data = res$data_long |>
                dplyr::filter(Gene %in% top_genes) |>
                dplyr::group_by(Gene) |>
                dplyr::slice_sample(n = 1) |>
                dplyr::ungroup(),
            aes(label = Gene, colour = Gene),
            size = 4, segment.size = 0.2, segment.color = "grey",
            max.overlaps = 50
        ) +
        scale_colour_manual(values = cols) +
        labs(x = "Ranked groups by expression", y = "Cumulative expression") +
        theme_minimal() +
        ZemmourLib::NoLegend() +
        ZemmourLib::NoGrid()
}

#' Skyline bar plot of one gene across groups
#'
#' The bar-per-cluster view used throughout Figures 7-8 and reproduced by the
#' public immgenT Skyline viewer.
#'
#' @param pseudobulk Seurat object of pseudobulk profiles.
#' @param gene Gene to plot.
#' @param group_levels Group order along the x axis.
#' @param group_colors Named colour vector indexed by `fill_by` values.
#' @param fill_by Optional named vector mapping group -> fill category (for
#'   example cluster -> lineage). When NULL each group is filled by its own name.
plot_gene_barplot <- function(pseudobulk, gene, group_levels,
                              group_colors, fill_by = NULL) {
    expr <- pseudobulk[["RNA"]]$data[gene, ]
    data <- data.frame(Group = names(expr), Expression = as.numeric(expr))
    data$Group <- factor(data$Group, levels = group_levels)
    data <- data[!is.na(data$Group), ]
    data$Fill <- if (is.null(fill_by)) {
        data$Group
    } else {
        factor(fill_by[as.character(data$Group)], levels = unique(fill_by))
    }

    ggplot(data, aes(x = Group, y = Expression, fill = Fill)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = group_colors) +
        labs(title = gene, x = NULL, y = "Expression log1p CP10K") +
        theme_minimal() +
        theme(
            axis.text.x  = element_text(size = 8, angle = 90,
                                        hjust = 1, vjust = 0.5),
            axis.title.y = element_text(size = 10)
        ) +
        ZemmourLib::NoGrid() +
        ZemmourLib::NoLegend()
}

#' Column annotation table and palettes for the pseudobulk heatmaps
#'
#' @param pseudobulk Seurat object of pseudobulk profiles.
#' @param cols Column names of the matrix being plotted.
#' @param ann_cols Metadata columns to annotate with.
#' @param ann_palettes Named list of colour vectors, one per annotation column.
#' @return List with `annotation_col` and `annotation_colors` for pheatmap.
heatmap_annotation <- function(pseudobulk, cols, ann_cols, ann_palettes = NULL) {
    if (is.null(ann_cols) || !length(ann_cols)) {
        return(list(annotation_col = NULL, annotation_colors = NULL))
    }
    md <- pseudobulk@meta.data
    if (!all(cols %in% rownames(md))) {
        stop("Could not match pseudobulk columns to pseudobulk@meta.data rows.")
    }
    md2 <- md[cols, , drop = FALSE]
    ann_cols_use <- intersect(ann_cols, colnames(md2))
    list(
        annotation_col = md2[, ann_cols_use, drop = FALSE],
        annotation_colors = if (is.null(ann_palettes)) NULL else
            ann_palettes[names(ann_palettes) %in% ann_cols_use]
    )
}

#' Expression heatmap of the most specific genes
#'
#' @param res Output of `gene_specificity_auc()`.
#' @param pseudobulk Seurat object of pseudobulk profiles.
#' @param n_top Number of top-scoring genes to show.
#' @param col_order Optional fixed column order; when supplied, columns are not
#'   clustered.
#' @param colours Heatmap colour ramp.
#' @param ann_cols,ann_palettes Passed to `heatmap_annotation()`.
#' @return A pheatmap object.
plot_gene_heatmap <- function(res, pseudobulk, n_top = 20,
                              col_order = NULL,
                              colours = ColorRamp,
                              ann_cols = c("annotation_level1",
                                           "annotation_level2_group"),
                              ann_palettes = NULL,
                              main = "Expression log1p CP10K",
                              cluster_rows = TRUE) {
    top_genes <- head(res$auc$Gene, n_top)
    mat <- res$mat[top_genes, , drop = FALSE]
    if (!is.null(col_order)) {
        mat <- mat[, intersect(col_order, colnames(mat)), drop = FALSE]
    }
    ann <- heatmap_annotation(pseudobulk, colnames(mat), ann_cols, ann_palettes)

    pheatmap::pheatmap(
        mat,
        color = colours,
        cluster_rows = cluster_rows,
        cluster_cols = is.null(col_order),
        main = main,
        annotation_col = ann$annotation_col,
        annotation_colors = ann$annotation_colors,
        silent = TRUE
    )
}
