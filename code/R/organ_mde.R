# Baseline all-T MDE gallery, one plot per organ.
#
# Shared by Figure 6a and Extended Data Figure 8, which are the same plot drawn
# for two disjoint sets of organs: Figure 6a covers the twelve well-sampled
# sites and Extended Data Figure 8 the five sparsely sampled ones.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})

# Organs shown by Figure 6a, in panel order.
ORGANS_FIGURE6 <- c(
    "blood", "spleen", "LN",
    "colon LP", "colon epi", "small intestine LP", "small intestine epi",
    "lung", "kidney", "liver", "skin", "mammary gland"
)

# Organs shown by Extended Data Figure 8, in panel order.
ORGANS_FIGURES8 <- c(
    "submandibular gland", "uterus", "peritoneal cavity", "placenta",
    "bone marrow"
)

# CNS is deliberately in neither list: at baseline it contributes too few cells
# for a per-organ embedding to be read, and it appears in the paper only
# through the composition dot plots and the alluvial plot.

#' Baseline cells, for the per-organ MDE gallery
#'
#' organ_simplified0 rather than organ_simplified: the latter pools sites that
#' the gallery shows separately.
#'
#' @param so_orig The full atlas object.
#' @return Healthy, non-thymic cells.
baseline_cells <- function(so_orig) {
    so_orig[, so_orig@meta.data %>%
                filter(condition_broad == "healthy",
                       organ_simplified0 != "thymus") %>%
                rownames()]
}

#' Write one organ's baseline MDE panel
#'
#' Organs are capped at 10,000 cells so the visual density of a panel reflects
#' composition rather than how deeply that tissue was sequenced. Callers should
#' set a seed.
#'
#' @param so_baseline Output of `baseline_cells()`.
#' @param organ One value of `organ_simplified0`.
#' @param figure_dir Figure directory, e.g. "Figure 6".
#' @param prefix Panel file name prefix, e.g. "6a_MDE_baseline".
plot_organ_mde <- function(so_baseline, organ, figure_dir, prefix,
                           max_cells = 10000) {
    cells <- colnames(so_baseline)[so_baseline$organ_simplified0 %in% organ]
    if (!length(cells)) {
        warning("no baseline cells for organ '", organ, "'; skipping")
        return(invisible(NULL))
    }
    if (length(cells) > max_cells) {
        cells <- sample(cells, size = max_cells, replace = FALSE)
    }
    panel_pdf(figure_dir,
              sprintf("%s_%s", prefix, gsub("[ /]", "", organ)), 6, 6)
    print(
        DimPlot(object = so_baseline[, cells],
                reduction = "mde2_totalvi_20241006",
                pt.size = 0.5, raster = FALSE,
                group.by = "annotation_level1", cols = mypal_level1,
                alpha = 1) +
            ggtitle(sprintf("%s (n = %s)", organ,
                            format(length(cells), big.mark = ","))) +
            ZemmourLib::NoLegend()
    )
    dev.off()
}
