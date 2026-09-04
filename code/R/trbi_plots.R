# Plotting helpers for the T-RBI integration panels (Figure 5, Extended Data
# Figure 6).
#
# Every T-RBI panel is the same two-layer construction: the whole immgenT atlas
# drawn in grey as a positional reference, with the mapped query cells drawn on
# top. The atlas layer is built from the atlas MDE embedding and the query layer
# from the anchored MDE coordinates T-RBI returns, so the two are on identical
# axes by construction. Both layers use scattermore, which rasterises the point
# cloud; at ~680,000 background points a vector layer produces PDFs too large
# to open.

suppressPackageStartupMessages({
    library(ggplot2)
    library(scattermore)
    library(dplyr)
})

#' Grey atlas background for a T-RBI panel
#'
#' @param so_atlas The immgenT atlas Seurat object.
#' @param reduction Atlas embedding to draw.
#' @return A ggplot with a single grey scattermore layer, to be added to.
trbi_background <- function(so_atlas, reduction = "mde2_totalvi_20241006") {
    emb <- so_atlas[[reduction]]@cell.embeddings
    df <- data.frame(dim1 = emb[, 1], dim2 = emb[, 2])
    ggplot(df) +
        geom_scattermore(aes(dim1, dim2), colour = "grey", alpha = 0.5,
                         pixels = c(1024, 1024))
}

#' Query-cell coordinates and metadata as a plain data frame
#'
#' @param so_query Merged Seurat object of the mapped external studies.
#' @param reduction Anchored MDE reduction holding the query coordinates.
#' @return Data frame of the query metadata with `dim1` / `dim2` appended.
trbi_foreground_df <- function(so_query, reduction = "mde_incremental_allT") {
    emb <- so_query[[reduction]]@cell.embeddings
    data.frame(
        so_query@meta.data[rownames(emb), ],
        dim1 = emb[, 1],
        dim2 = emb[, 2]
    )
}

# The atlas MDE occupies roughly [-2.5, 2.5] on both axes. Fixing the limits
# keeps every T-RBI panel on the same scale, so galleries can be compared
# panel to panel and against the atlas figures.
TRBI_XLIM <- c(-2.5, 2.5)
TRBI_YLIM <- c(-2.5, 2.5)

#' Standard theme for the T-RBI MDE panels
trbi_theme <- function() {
    list(
        theme_void(),
        xlim(TRBI_XLIM),
        ylim(TRBI_YLIM),
        ZemmourLib::NoGrid()
    )
}
