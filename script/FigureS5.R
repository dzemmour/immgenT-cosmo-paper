# Extended Data Figure 5. Distribution of stress-related transcriptional
# signatures across the immgenT atlas.
#
# Panels produced:
#   S5a-S5h  Single-cell scores for eight representative signatures on the
#            all-T MDE.
#   S5i      Mean signature scores across immune perturbations.
#   S5j      Mean signature scores across anatomical sites.
#   S5k      Mean signature scores across clusters.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Stress signatures" (line 2965). The Rmd
# read signatures_stress.Rds and signature_stress_scores.csv from out_dir rather
# than data_dir, which only worked because both directories held copies; both
# are inputs and are read from data/ here.
#
# The Rmd's FeaturePlot loop drew all fourteen signatures; the eight published
# panels are selected explicitly below, in panel order.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   signatures_stress.Rds               [curated input]
#   signature_stress_scores.csv         [code/pipeline/02]

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ComplexHeatmap)
    library(circlize)
})

source("code/R/setup.R")

figure_dir <- "Extended Data Figure 5"

# The eight signatures shown as MDE panels, in published order.
SIGNATURE_PANELS <- c(
    "S5a" = "dissociation1h",
    "S5b" = "dissociation2h",
    "S5c" = "ier",
    "S5d" = "isr",
    "S5e" = "erstress_GO0034976",
    "S5f" = "isg_immgen",
    "S5g" = "hypoxia_GO0001666",
    "S5h" = "mt_genes"
)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
sig_list <- readRDS(sprintf("%s/signatures_stress.Rds", data_path))

# Scores were computed once with AddModuleScore() on the complete object and
# cached, because scoring fourteen signatures against 683,000 cells with 100
# control genes each is the expensive step. See code/pipeline/02_stress_scores.R.
sig_scores <- read.csv(sprintf("%s/signature_stress_scores.csv", data_path),
                       stringsAsFactors = FALSE, check.names = FALSE)
sig_scores <- sig_scores[match(colnames(so_orig), sig_scores$cellID), ]
stopifnot(identical(sig_scores$cellID, colnames(so_orig)))

score_cols <- setdiff(colnames(sig_scores), "cellID")
so_orig@meta.data[, score_cols] <- sig_scores[, score_cols, drop = FALSE]

so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus") %>% rownames()]

# ============================================================
# S5a-S5h: Signature scores on the all-T MDE
# ============================================================
# Capped at the 5th and 99th percentiles: module scores have long tails in both
# directions, and without the cap a few extreme cells set the scale and
# everything else renders as one flat colour.
for (panel in names(SIGNATURE_PANELS)) {
    sig <- SIGNATURE_PANELS[[panel]]
    panel_pdf(figure_dir, sprintf("%s_MDE_%s", panel, sig), 5, 5)
    print(
        FeaturePlot(
            so, reduction = "mde2_totalvi_20241006", features = sig,
            order = TRUE, raster = TRUE, raster.dpi = c(512, 512),
            min.cutoff = "q5", max.cutoff = "q99"
        ) +
            theme_void() + ggtitle(sig) +
            xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + NoLegend()
    )
    dev.off()
}

# ============================================================
# S5i-S5k: Mean signature scores by condition, organ and cluster
# ============================================================
# Rows are mean-centred, so the colour is a signature's deviation from its own
# average across the groupings shown. Without centring the plot would be
# dominated by the fact that some signatures have larger scores than others
# everywhere, which says nothing about where stress is elevated.
#
# ifn_gp16 and ifn_gp182 are excluded from the heatmaps: they are small
# interferon gene-programme subsets superseded by isg_immgen, which is shown.
heatmap_cols <- setdiff(names(sig_list), c("ifn_gp16", "ifn_gp182"))
heatmap_cols <- intersect(heatmap_cols, colnames(so@meta.data))

signature_heatmap <- function(meta, group_col, group_levels = NULL) {
    avg <- meta %>%
        filter(!is.na(.data[[group_col]])) %>%
        group_by(.data[[group_col]]) %>%
        summarise(across(all_of(heatmap_cols), ~ mean(.x, na.rm = TRUE)),
                  .groups = "drop")

    mat <- avg %>%
        column_to_rownames(group_col) %>%
        as.matrix() %>%
        t()

    if (!is.null(group_levels)) {
        mat <- mat[, intersect(group_levels, colnames(mat)), drop = FALSE]
    }
    mat <- sweep(mat, 1, rowMeans(mat, na.rm = TRUE), "-")

    lim <- max(abs(mat), na.rm = TRUE)
    Heatmap(
        mat,
        name = "Mean score (centered)",
        col = colorRamp2(c(-lim, 0, lim), c("#3B4CC0", "white", "#B40426")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_rot = 90
    )
}

# By perturbation. Computed on the full object, including thymus, because the
# thymic samples are their own baseline condition. "other" pools perturbations
# too sparse to interpret.
meta_cond <- so_orig@meta.data %>%
    filter(condition_detailed_simplified != "other")
panel_pdf(figure_dir, "S5i_heatmap_by_condition", 15, 5)
draw(signature_heatmap(meta_cond, "condition_detailed_simplified"))
dev.off()

# By anatomical site. SLO is dropped: it is a catch-all for secondary lymphoid
# tissue that does not resolve to a specific organ.
meta_organ <- so@meta.data %>% filter(organ_simplified != "SLO")
panel_pdf(figure_dir, "S5j_heatmap_by_organ", 6, 4)
draw(signature_heatmap(meta_organ, "organ_simplified",
                       levels(so$organ_simplified)))
dev.off()

# By cluster, in lineage order.
panel_pdf(figure_dir, "S5k_heatmap_by_cluster", 18, 3)
draw(signature_heatmap(so@meta.data, "annotation_level2",
                       levels(so$annotation_level2)))
dev.off()
