# Discovery score for query cells T-RBI leaves unclassified.
#
# The score answers whether unclassified query cells form a coherent population
# of their own -- a candidate state missing from the reference -- or are
# scattered among confidently annotated cells at cluster boundaries. It is
# deliberately computed in the query dataset's own PCA space, independently of
# the immgenT integration, so that it cannot inherit artefacts of the mapping it
# is meant to audit.
#
# These functions produced the cached score table read by script/Figure5.R and
# script/FigureS6.R; they are retained here to document its provenance.

suppressPackageStartupMessages({
    library(RANN)
    library(dplyr)
})

#' Nearest-neighbour novelty score for query cells
#'
#' For each unclassified cell, compares its mean distance to the k nearest
#' annotated ("reference") cells with its mean distance to the k nearest other
#' unclassified ("query") cells. A ratio above 1 means the cell sits closer to
#' other unclassified cells than to annotated ones.
#'
#' @param old_mat Embedding of the annotated cells.
#' @param new_mat Embedding of the unclassified cells, same dimensions.
#' @param k Number of neighbours.
#' @return Tibble with `mean_d_ref`, `mean_d_qry`, `ratio` and two rescaled
#'   forms of the same contrast (`score_log`, `score_sym`).
knn_novelty_scores <- function(old_mat, new_mat, k = 10) {
    old_mat <- as.matrix(old_mat)
    new_mat <- as.matrix(new_mat)
    stopifnot(ncol(old_mat) >= 2, ncol(new_mat) >= 2,
              ncol(old_mat) == ncol(new_mat))
    if (nrow(new_mat) <= k) stop("Need k < number of new points.")

    ref_knn <- nn2(data = old_mat, query = new_mat, k = k)
    mean_d_ref <- rowMeans(ref_knn$nn.dists)

    # k + 1 neighbours, then drop the first column: querying a point against
    # its own set returns itself at distance 0.
    qry_knn <- nn2(data = new_mat, query = new_mat, k = k + 1)
    mean_d_qry <- rowMeans(qry_knn$nn.dists[, -1, drop = FALSE])

    ratio <- mean_d_ref / mean_d_qry

    tibble::tibble(
        idx = seq_len(nrow(new_mat)),
        mean_d_ref = mean_d_ref,
        mean_d_qry = mean_d_qry,
        ratio = ratio,
        score_log = log(ratio),
        score_sym = (mean_d_ref - mean_d_qry) / (mean_d_ref + mean_d_qry)
    )
}

#' Summarise a novelty score table
#'
#' @param scores_tbl Output of `knn_novelty_scores()`.
#' @param ratio_threshold Ratio above which a cell is counted as a candidate
#'   novel state. The paper uses a conservative 1.1.
summarize_discovery <- function(scores_tbl, ratio_threshold = 2) {
    tibble::tibble(
        k = NA_integer_,
        ratio_threshold = ratio_threshold,
        discovery_ncells = sum(na.omit(scores_tbl$ratio > ratio_threshold)),
        median_ratio = median(na.omit(scores_tbl$ratio)),
        q90_ratio = quantile(na.omit(scores_tbl$ratio), 0.9),
        mean_log_score = mean(na.omit(scores_tbl$score_log))
    )
}
