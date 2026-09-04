# Shared data loading and colour palettes.
#
# Every script in script/ begins by sourcing this file. It defines the colour
# palettes used throughout the paper and the two loaders that return the atlas
# object in the two forms the figure scripts need.
#
# Palettes come from the ZemmourLib package (immgent_colors) so that a cluster
# or lineage is drawn in the same colour on every panel of every figure. The
# per-experiment (IGT) palette is the one exception: it is a fixed assignment
# stored in data/color_palette_igt.csv, so that experiment IGT5 keeps its colour
# even as experiments are added to or removed from a plot.

suppressPackageStartupMessages({
    library(Seurat)
    library(ZemmourLib)
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(pals)
    library(scales)
})

data_path <- "data"

# --- Colour palettes ---------------------------------------------------------

# Diverging ramp used for the pseudobulk expression heatmaps.
ColorRamp <- rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100))

# Cluster (level 2), lineage (level 1) and organ palettes. "not classified" is
# added to the first two because T-RBI query cells that fail the confidence
# threshold are plotted alongside annotated cells (Figure 5, Extended Data
# Figure 6) and must be visually distinct from every real cluster.
mypal_level2 <- ZemmourLib::immgent_colors$level2
mypal_level2["not classified"] <- "black"
mypal_level1 <- ZemmourLib::immgent_colors$level1
mypal_level1["not classified"] <- "black"
mypal_organ <- ZemmourLib::immgent_colors$organ_simplified

# Activation state annotation used on the effector and transcription factor
# heatmaps.
mypal_level2group <- c(
    "resting"       = "blue",
    "activated"     = "red",
    "miniverse"     = "darkgreen",
    "proliferating" = "black",
    "other"         = "grey",
    "preT"          = "grey"
)

#' Fixed per-experiment (IGT) colour palette
#'
#' @return Named character vector of colours, one per IGT identifier.
load_igt_palette <- function(data_dir = data_path) {
    tmp <- read.table(
        sprintf("%s/color_palette_igt.csv", data_dir),
        header = FALSE, sep = ","
    )
    setNames(tmp[, 2], tmp[, 1])
}

# Large categorical palette for variables with more levels than any curated
# palette covers (experiment identifiers, external study names, per-experiment
# transcriptomic clusters). Glasbey and Polychrome first because they are
# designed to maximise perceptual distance, then the ColorBrewer qualitative
# palettes concatenated. The fourth ColorBrewer colour is dropped because it is
# a near-white that disappears against the panel background.
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
mypal1 <- unique(unlist(mapply(
    brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)
)))
mypal1 <- mypal1[-4]
mypal <- c(glasbey(), polychrome(), mypal1)
names(mypal) <- NULL

# --- Atlas loaders -----------------------------------------------------------

#' Load the immgenT atlas
#'
#' Reads `immgenT_seurat_ADT_GeneSubset.Rds`, the gene- and protein-subsetted
#' atlas object (see code/README.md for how it is derived from the complete
#' object). RNA counts in this object are **already log-normalised**; do not
#' call `NormalizeData()` on the RNA assay again. The ADT assay is not
#' normalised, and each script that plots protein expression normalises it
#' locally.
#'
#' @param drop_thymus Exclude thymic cells. Every published panel except the
#'   double-positive lineage MDE is restricted to non-thymic cells, because
#'   thymocyte states are the subject of a companion manuscript.
#' @param citeseq_only Restrict to cells with CITE-seq data, required for any
#'   panel showing surface protein.
#' @return A Seurat object.
load_immgent <- function(data_dir = data_path,
                         drop_thymus = TRUE,
                         citeseq_only = FALSE) {
    so <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_dir))
    keep <- so@meta.data
    if (drop_thymus) keep <- keep %>% filter(organ_simplified != "thymus")
    if (citeseq_only) keep <- keep %>% filter(cite_seq == TRUE)
    so[, rownames(keep)]
}

#' A single panel from MyDimPlotHighlight
#'
#' `ZemmourLib::MyDimPlotHighlight()` both prints its plots (governed by
#' `print_plot1` / `print_plot2`) and returns all of them, so calling it and
#' printing the result writes two or three pages to the open device. This
#' suppresses the internal printing and hands back only the variant asked for:
#'
#'   plot1  titled, with axes and legend (theme_minimal)
#'   plot2  bare: theme_void, no legend, no axes
#'   plot3  plot2 plus a title and the cluster labels placed at each cluster's
#'          median position; NULL unless `labelclusters = TRUE`
#'
#' @param which One of "plot1", "plot2", "plot3".
#' @param ... Passed to `MyDimPlotHighlight()`.
#' @return A ggplot.
highlight_mde <- function(..., which = c("plot2", "plot1", "plot3")) {
    which <- match.arg(which)
    p <- ZemmourLib::MyDimPlotHighlight(
        ..., print_plot1 = FALSE, print_plot2 = FALSE
    )
    if (is.null(p[[which]])) {
        stop("MyDimPlotHighlight() returned no ", which,
             " (plot3 requires labelclusters = TRUE)")
    }
    p[[which]]
}

#' Open a PDF device for one panel
#'
#' Panels are written one PDF per panel into `figures/<figure>/<panel>.pdf`.
#' `useDingbats = FALSE` keeps point glyphs as vector paths rather than
#' Dingbats characters, which some vector editors fail to render.
#'
#' @param figure Figure directory name, e.g. "Figure 4".
#' @param panel Panel file name without extension, e.g. "4e".
panel_pdf <- function(figure, panel, width, height, dir = "figures") {
    outdir <- file.path(dir, figure)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(
        file.path(outdir, paste0(panel, ".pdf")),
        width = width, height = height, useDingbats = FALSE
    )
}
