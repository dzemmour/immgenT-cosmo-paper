# Figure 10. AI-based reference mapping reveals the limits of marker- and
# signature-defined T-cell identity.
#
# Panels produced:
#   10a  Trm cells (P14, 30+ dpi LCMV Armstrong, gut) on the CD8-specific MDE.
#   10b  CD103+CD69+ cells across the whole CD8 MDE.
#   10c  CD69 vs CD103 protein in memory cells from small intestine, prostate,
#        salivary gland and spleen, coloured by cluster.
#   10f  Alluvial plot of CD8.Q against organ and immune perturbation.
#
# Panels 10d (a published Trm gene signature scored across the atlas) and 10e
# (T-RBI mapping of the cells that signature was derived from) both require the
# Milner et al. signature and dataset, neither of which is among the cached
# inputs; see analysis/Figure10.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 9" (line 2442), whose numbering is
# from an earlier draft; the Trm material is Figure 10 in the published paper.
# The Rmd's chunks under that header also produce Extended Data Figure 11c (the
# Treg core signature), which lives in script/FigureS11.R.
#
# Panel 10a: the Rmd drew two variants of this plot, one over four non-lymphoid
# tissues and one over the gut alone. The published panel is the gut-only one.
#
# Panel 10b: the published panel renders the gated cells as a white-to-purple
# cell-density map with a Low-High legend, whereas the Rmd marks them as solid
# blue points. The Rmd's version is kept: the two select the same cells and
# agree panel-for-panel on where they sit, and MyDimPlotHighlightDensity()'s
# densCols() colouring is dominated by the tight CD8.Q ball, which pushes the
# broadly distributed cells -- the panel's whole point -- to the white end of
# any ramp. The exact ramp behind the published rendering is not recorded.
#
# Panel 10c: the Rmd drew six CD69/CD103 scatters -- the four published here
# plus one pooling all four non-lymphoid tissues and one coloured by organ
# rather than by cluster. The latter two are not published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggalluvial)
    library(tidyr)
    library(viridis)
})

source("code/R/setup.R")
source("code/R/alluvial.R")

figure_dir <- "Figure 10"
set.seed(1)

# The two experiments that profiled P14 (LCMV GP33-specific) CD8 cells at
# memory timepoints across tissues.
TRM_EXPERIMENTS <- c("IGT38", "IGT40")

# Protein positivity thresholds for the two tissue-residency markers, read off
# the bimodal distributions in panel c.
CD103_GATE <- 4.5
CD69_GATE <- 4
# Panel b gates CD103 one step lower than panel c's dashed line, as the
# original analysis did: the panel is making the case that the marker
# combination is not specific, so the more permissive gate is the harder test.
CD103_GATE_DENSITY <- 4

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
cd8 <- so_orig[, so_orig$annotation_level1 == "CD8"]

# ============================================================
# 10a: Trm cells on the CD8-specific MDE
# ============================================================
# Bona fide memory cells, identified by antigen specificity and timepoint
# rather than by marker expression, so that where they land is an independent
# check on the reference rather than a restatement of the gate.
trm_cells <- cd8@meta.data %>%
    filter(Ag_spe == "P14",
           IGT %in% TRM_EXPERIMENTS,
           organ_simplified %in% c("small intestine epi",
                                   "small intestine LP")) %>%
    rownames()

panel_pdf(figure_dir, "10a_CD8MDE_P14_Trm_gut", 5, 5)
print(highlight_mde(
    seurat_object = cd8, umap_to_plot = "mde_incremental",
    cells_to_highlight = trm_cells, background_alpha = 0.5,
    highlight_column_name = "annotation_level2", pixels = c(512, 512),
    mycols = mypal_level2,
    title = "Trm cells (IGT38, 40)\nP14 cells, 30+ dpi LCMVarm, gut",
    highlight_size = 2, highlight_alpha = 0.5, labelclusters = FALSE,
    which = "plot1"
))
dev.off()

# ============================================================
# 10b: CD103+CD69+ cells across the CD8 MDE
# ============================================================
# The complement of panel a: instead of asking where known memory cells land,
# this asks where the canonical marker combination is found. It is not confined
# to the region panel a picks out, which is the figure's point.
cd8_cite <- cd8[, cd8$cite_seq == TRUE]
cd8_cite <- NormalizeData(cd8_cite, assay = "ADT",
                          normalization.method = "LogNormalize",
                          verbose = FALSE)

adt <- cd8_cite[["ADT"]]$data
cd8_cite$is_CD103p_CD69p <-
    adt["CD69", ] > CD69_GATE & adt["CD103", ] > CD103_GATE_DENSITY

panel_pdf(figure_dir, "10b_CD8MDE_CD103pos_CD69pos", 5, 5)
print(highlight_mde(
    seurat_object = cd8_cite, umap_to_plot = "mde_incremental",
    cells_to_highlight = cd8_cite@meta.data %>%
        filter(is_CD103p_CD69p) %>% rownames(),
    highlight_column_name = "is_CD103p_CD69p", pixels = c(512, 512),
    mycols = c("grey", "#0C03FF"), background_alpha = 0.5,
    highlight_size = 1, highlight_raster = TRUE, labelclusters = FALSE,
    which = "plot2"
))
dev.off()

# ============================================================
# 10c: CD69 vs CD103 by tissue, coloured by cluster
# ============================================================
# One panel per tissue rather than one pooled panel: the markers are variably
# expressed between tissues, and pooling averages that variation away. Cells
# are coloured by their immgenT cluster, so within a panel one can see that
# CD8.Q cells sit on both sides of the gate.
cd8_cite_all <- NormalizeData(cd8, assay = "ADT",
                              normalization.method = "LogNormalize",
                              verbose = FALSE)

tissue_panels <- list(
    "small_intestine" = list(
        organs = c("small intestine epi", "small intestine LP"),
        title = "Small intestine"
    ),
    "prostate" = list(organs = "prostate", title = "Prostate"),
    "salivary_gland" = list(organs = "submandibular gland",
                            title = "Salivary gland"),
    "spleen" = list(organs = c("spleen", "LN"), title = "Spleen")
)

for (nm in names(tissue_panels)) {
    spec <- tissue_panels[[nm]]
    cells <- cd8_cite_all@meta.data %>%
        filter(Ag_spe == "P14",
               IGT %in% TRM_EXPERIMENTS,
               organ_simplified %in% spec$organs) %>%
        rownames()

    p <- FeatureScatter(
        object = cd8_cite_all[, cells],
        feature1 = "CD103", feature2 = "CD69",
        group.by = "annotation_level2", cols = mypal_level2,
        plot.cor = FALSE, pt.size = 2
    ) +
        ggtitle(spec$title) +
        xlim(2, 7) + ylim(2, 6) +
        geom_hline(yintercept = CD69_GATE, linetype = "dashed",
                   colour = "grey") +
        geom_vline(xintercept = CD103_GATE, linetype = "dashed",
                   colour = "grey")

    panel_pdf(figure_dir, sprintf("10c_CD103_CD69_%s", nm), 5, 5)
    print(p)
    dev.off()
}

rm(cd8_cite, cd8_cite_all)

# ============================================================
# 10f: CD8.Q against organ and immune perturbation
# ============================================================
# Perturbations are downsampled to 1,000 cells each so that flow width
# reflects how broadly a cluster is distributed rather than how many cells a
# given experiment contributed. CD8.P and CD8.wM (the proliferating and
# miniverse clusters) and SLO are excluded: they draw cells from every
# condition and would obscure the rest.
df0 <- so_orig@meta.data %>%
    filter(annotation_level1 == "CD8",
           !(annotation_level2 %in% c("CD8.P", "CD8.wM")),
           !(organ_simplified %in% "SLO")) %>%
    group_by(condition_detailed_simplified) %>%
    slice_sample(n = 1000) %>%
    ungroup()

# Each axis is ordered by hierarchical clustering of its own cross-tabulation,
# so that similar clusters, tissues and conditions sit next to each other and
# the ribbons cross as little as possible.
hclust_order <- function(df, rows, cols) {
    mat <- df %>%
        count(.data[[rows]], .data[[cols]], name = "count") %>%
        tidyr::pivot_wider(names_from = all_of(cols), values_from = count,
                           values_fill = 0) %>%
        as.data.frame()
    rownames(mat) <- mat[[rows]]
    mat[[rows]] <- NULL
    rownames(mat)[hclust(dist(mat))$order]
}

cluster_order <- hclust_order(df0, "annotation_level2", "organ_simplified")
organ_order <- hclust_order(df0, "organ_simplified", "annotation_level2")
cond_order <- hclust_order(df0, "condition_detailed_simplified",
                           "annotation_level2")

df3 <- df0 %>%
    count(organ_simplified, annotation_level2, condition_detailed_simplified,
          name = "count") %>%
    mutate(
        annotation_level2 = factor(annotation_level2, levels = cluster_order),
        organ_simplified = factor(organ_simplified, levels = organ_order),
        condition_detailed_simplified =
            factor(condition_detailed_simplified, levels = cond_order)
    )

highlight <- "CD8.Q"
df3_plot <- df3 %>%
    mutate(
        is_highlight = annotation_level2 %in% highlight,
        fill_key = ifelse(is_highlight, as.character(annotation_level2),
                          "other"),
        alpha_val = ifelse(is_highlight, 1, 0.25)
    )

panel_pdf(figure_dir, "10f_alluvial_CD8Q", 15, 5)
print(PlotAlluvialHighlight(
    df3_plot,
    mypal_alluv = c(other = "grey80", setNames("red", highlight)),
    title = sprintf("%s: organ -> cluster -> condition", highlight),
    alpha_range = c(0.2, 1)
))
dev.off()
