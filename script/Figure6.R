# Figure 6. Reuse and redistribution of shared T-cell states across tissues and
# immune perturbations.
#
# Panels produced:
#   6a  All-T MDE at baseline, one plot per organ, coloured by lineage.
#   6b  Cluster proportions across organs at baseline.
#   6c  Cluster proportions across immune perturbations in the colon.
#   6d  Cluster proportions across immune perturbations in the lung.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Fig 6a - S6 - MDE plots organ (Fig S6)"
# (line 1663) and "Figure 6 - Proportion of level2 ... in each organ"
# (line 1686). That first header's "S6" is stale: the supplementary counterpart
# of panel 6a is Extended Data Figure 8, and it is drawn by the same loop, so
# code/R/organ_mde.R holds the shared plotting code and the two organ lists,
# and script/FigureS8.R draws the five organs that figure shows. The Rmd looped
# over every baseline organ, which also emitted a CNS panel that neither figure
# publishes.
#
# The dot-plot section also produced an all-conditions version of 6b and a
# by-condition_broad version, neither of which is published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds                         [primary input]
#   annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds [primary input]
#   Sample_metadata_david_20260107_v11.csv                    [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source("code/R/setup.R")
source("code/R/composition.R")
source("code/R/organ_mde.R")

figure_dir <- "Figure 6"
set.seed(1)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus", cite_seq == TRUE) %>%
                  rownames()]

# ============================================================
# 6a: All-T MDE at baseline, per organ
# ============================================================
# Baseline only, so that what the panel shows is the resting composition of
# each tissue rather than the composition of whichever challenge happened to be
# applied there. The twelve organs and their order are fixed in
# code/R/organ_mde.R; the five sparsely sampled sites are Extended Data
# Figure 8.
so_baseline <- baseline_cells(so_orig)

for (organ in ORGANS_FIGURE6) {
    plot_organ_mde(so_baseline, organ, figure_dir, "6a_MDE_baseline")
}

rm(so_baseline)

# ============================================================
# 6b: Cluster proportions across organs at baseline
# ============================================================
# The proportion is within lineage, not within sample: it answers "of this
# tissue's CD8 cells, what fraction are CD8.Q", which is comparable between
# tissues that differ wildly in lineage composition. Restricted to samples
# sorted on all T cells or all CD45+ cells, because a sample sorted on, say,
# Tregs cannot report a lineage's internal composition.
m3 <- readRDS(sprintf(
    "%s/annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds", data_path
))

m_healthy <- m3 %>%
    filter(organ_simplified != "thymus",
           target_cells_simplified %in% c("CD45p", "allT"),
           condition_broad == "healthy") %>%
    select(-contains("thymocyte"))

res <- make_comp_dotdata_and_plot(data = m_healthy, so = so,
                                  condition_col = "organ_simplified")

panel_pdf(figure_dir, "6b_DotPlot_level2_by_organ_baseline", 25, 5)
print(
    res$plot +
        scale_fill_gradientn(
            colours = colorRampPalette(c("white", "black"))(10),
            name = "Mean prop"
        )
)
dev.off()

# ============================================================
# 6c: Cluster proportions across perturbations in the colon
# ============================================================
# "other" pools perturbations too sparsely sampled to interpret, and
# SFB_pregnancy confounds two perturbations at once; both are excluded rather
# than shown as a row that cannot be read.
m_colon <- m3 %>%
    filter(organ_simplified == "colon LP",
           target_cells_simplified %in% c("CD45p", "allT"),
           !(condition_detailed_simplified %in% c("other", "SFB_pregnancy"))) %>%
    select(-contains("thymocyte"))

res <- make_comp_dotdata_and_plot(
    data = m_colon, so = so,
    condition_col = "condition_detailed_simplified"
)

panel_pdf(figure_dir, "6c_DotPlot_level2_by_condition_colon", 22, 3.5)
print(
    res$plot +
        scale_fill_gradientn(
            colours = colorRampPalette(c("white", mypal_organ["colon"]))(10),
            name = "Mean prop"
        )
)
dev.off()

# ============================================================
# 6d: Cluster proportions across perturbations in the lung
# ============================================================
m_lung <- m3 %>%
    filter(organ_simplified == "lung",
           target_cells_simplified %in% c("CD45p", "allT"),
           condition_detailed_simplified != "other") %>%
    select(-contains("thymocyte"))

res <- make_comp_dotdata_and_plot(
    data = m_lung, so = so,
    condition_col = "condition_detailed_simplified"
)

panel_pdf(figure_dir, "6d_DotPlot_level2_by_condition_lung", 22, 3.5)
print(res$plot)
dev.off()
