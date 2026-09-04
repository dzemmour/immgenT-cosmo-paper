# Cluster robustness metrics: silhouette conservation across integration, and
# centroid-based cluster separability.
#
# Both quantities are computed once per experiment (IGT) and cached, because
# each requires an independent clustering or PCA of that experiment's cells.
# script/Figure4.R and script/FigureS4.R read the cached .Rds files; the
# functions here are what produced them, retained so those files have a
# documented provenance. See code/pipeline/01_cluster_metrics.R for the driver.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cluster)
})

#' Silhouette width per cluster, before and after integration
#'
#' Clusters an experiment's cells de novo on its own transcriptome, so the
#' clustering is independent of the atlas annotation, then scores how well
#' separated those clusters are in two spaces: the experiment's own PCA
#' ("pre-integration") and the shared totalVI latent space
#' ("post-integration").
#'
#' @param so Seurat object holding a single experiment's cells.
#' @param reduc_integrated Name of the integrated latent-space reduction.
#' @return One row per de novo cluster, with mean and median silhouette width
#'   pre- and post-integration.
BioConservationSilhouette <- function(so, reduc_integrated) {
    message("Clustering within the experiment...")
    so <- so %>%
        NormalizeData(normalization.method = "LogNormalize",
                      scale.factor = 1e4, verbose = FALSE) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000,
                             verbose = FALSE) %>%
        ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
        RunPCA(features = VariableFeatures(.), npcs = 50, verbose = FALSE) %>%
        FindNeighbors(dims = 1:30, verbose = FALSE) %>%
        FindClusters(resolution = 1, verbose = FALSE)

    summarise_sil <- function(embeddings, clusters) {
        sil <- silhouette(clusters, dist = dist(embeddings))
        as.data.frame(sil) %>%
            group_by(cluster) %>%
            summarise(
                mean_sil_width      = mean(sil_width),
                median_sil_width    = median(sil_width),
                total               = n(),
                positive_count      = sum(sil_width > 0),
                percentage_positive = (positive_count / total) * 100,
                .groups = "drop"
            ) %>%
            as.data.frame()
    }

    cl <- as.integer(so[["RNA_snn_res.1"]][, 1])

    message("Calculating silhouette pre integration...")
    sil_summary_pre <- summarise_sil(so[["pca"]]@cell.embeddings, cl)

    message("Calculating silhouette post integration...")
    sil_summary_post <- summarise_sil(so[[reduc_integrated]]@cell.embeddings, cl)

    merge(sil_summary_pre, sil_summary_post,
          by = "cluster", suffixes = c(".pre", ".post"))
}

#' Centroid-based cluster separability
#'
#' For each cluster, compares the mean distance from its cells to its own
#' centroid, a(C), with the distance from that centroid to the other cluster
#' centroids, b(C), as s = (b - a) / max(a, b). The score is bounded in
#' [-1, 1]; positive values mean cells of a cluster are on average closer to
#' each other than to other clusters.
#'
#' Used in place of a per-cell silhouette because it is insensitive to cluster
#' size and to local neighbourhood structure.
#'
#' @param emb Cell x dimension embedding matrix.
#' @param clusters Cluster label per cell, in the row order of `emb`.
#' @param b_mode "mean" compares each centroid to all other centroids;
#'   "min" compares to the nearest other centroid only.
#' @return One row per cluster with columns `cluster`, `a`, `b`, `s`.
cluster_separability_centroid <- function(emb, clusters,
                                          b_mode = c("mean", "min")) {
    b_mode <- match.arg(b_mode)

    emb <- as.matrix(emb)
    if (nrow(emb) != length(clusters)) {
        stop("length(clusters) must match nrow(emb)")
    }
    clusters <- as.factor(clusters)
    levs <- levels(clusters)

    centroids <- vapply(levs, function(cl) {
        colMeans(emb[clusters == cl, , drop = FALSE])
    }, FUN.VALUE = numeric(ncol(emb)))
    centroids <- t(centroids)
    rownames(centroids) <- levs

    # a(C): mean distance of member cells to their own centroid.
    a <- vapply(levs, function(cl) {
        X <- emb[clusters == cl, , drop = FALSE]
        mu <- centroids[cl, , drop = FALSE]
        mean(sqrt(rowSums((t(apply(X, 1, function(r) r - mu)))^2)))
    }, FUN.VALUE = numeric(1))
    names(a) <- levs

    Dcc <- as.matrix(dist(centroids))

    # b(C): mean (or minimum) distance from this centroid to the others.
    b <- vapply(levs, function(cl) {
        others <- setdiff(levs, cl)
        if (length(others) == 0) return(NA_real_)
        if (b_mode == "mean") mean(Dcc[cl, others]) else min(Dcc[cl, others])
    }, FUN.VALUE = numeric(1))
    names(b) <- levs

    s <- (b - a) / pmax(a, b)
    data.frame(
        cluster = levs,
        a = unname(a),
        b = unname(b),
        s = unname(s),
        stringsAsFactors = FALSE
    )
}
