# Extended Data Figure 11. Markers and gene signatures across immgenT states.
#
# Panels produced:
#   S11a  CD4-specific MDE highlighting the Th1, Th2 and Th17 "tip" clusters.
#   S11b  Th-defining receptor transcripts on the CD4-specific MDE.
#   S11c  The core Treg gene signature scored across the whole atlas, on the
#         all-T MDE.
#
# --- internal ---
# Panel c is ported from the "Treg signature score from Zemmour et al. 2021"
# chunk of cosmo_paper.Rmd section "Figure 9" (line 2442). The score is read
# from cache rather than recomputed with AddModuleScore() as the Rmd does; the
# code that produced it is code/pipeline/02_signature_scores.R.
#
# Panels a and b have no counterpart in the Rmd. They are reconstructed here
# from the cluster assignments and the marker transcripts the legend names.
# They are the one place in this repository that needs the complete atlas
# object rather than the gene-subsetted one: three of the six transcripts
# (Cxcr3, Ccr4, Il4ra) were dropped by the DietSeurat step, and no panel in the
# Rmd needed them.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   immgenT_seurat_ADT_complete.Rds     [primary input]
#   treg_core_score.txt                 [code/pipeline/02]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source("code/R/setup.R")

figure_dir <- "Extended Data Figure 11"

# The polarised CD4 "tip" clusters, named in the figure legend.
TH_TIPS <- list(
    Th1  = c("CD4.R", "CD4.S"),
    Th2  = c("CD4.T", "CD4.U"),
    Th17 = c("CD4.X", "CD4.Y")
)

TH_MARKERS <- c("Cxcr3", "Il12rb2", "Ccr4", "Il4ra", "Ccr6", "Il23r")

# ============================================================
# S11c: Core Treg signature across the atlas
# ============================================================
# The cached score, over every T cell rather than only Tregs. Capped at the 5th
# and 90th percentiles for the same reason as the Extended Data Figure 5
# panels.
so <- load_immgent()

treg_core <- read.table(sprintf("%s/treg_core_score.txt", data_path),
                        header = TRUE, stringsAsFactors = FALSE)
so$treg_core <- treg_core[match(colnames(so), rownames(treg_core)), 1]

panel_pdf(figure_dir, "S11c_MDE_treg_core_signature", 5, 5)
print(
    FeaturePlot(so, reduction = "mde2_totalvi_20241006", features = "treg_core",
                order = TRUE, raster = TRUE, raster.dpi = c(512, 512),
                min.cutoff = "q5", max.cutoff = "q90") +
        ggtitle("Treg core signature") +
        xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
)
dev.off()

rm(so)
invisible(gc())

# ============================================================
# S11a: Th tip clusters on the CD4 MDE
# ============================================================
# Loaded from the complete object because panel b needs transcripts the
# gene-subsetted object does not carry; panel a uses the same object so both
# sit on identical axes.
so_full <- readRDS(sprintf("%s/immgenT_seurat_ADT_complete.Rds", data_path))
cd4 <- so_full[, so_full@meta.data %>%
                   filter(annotation_level1 == "CD4",
                          organ_simplified != "thymus") %>%
                   rownames()]
rm(so_full)
invisible(gc())

cd4$th_tip <- NA_character_
for (tip in names(TH_TIPS)) {
    cd4$th_tip[cd4$annotation_level2 %in% TH_TIPS[[tip]]] <- tip
}

panel_pdf(figure_dir, "S11a_CD4MDE_Th_tips", 5, 5)
print(highlight_mde(
    seurat_object = cd4, umap_to_plot = "mde_incremental",
    cells_to_highlight = colnames(cd4)[!is.na(cd4$th_tip)],
    highlight_column_name = "th_tip", pixels = c(512, 512),
    mycols = c(Th1 = "darkorange2", Th2 = "deeppink", Th17 = "chartreuse3"),
    title = "Th1 (CD4.R/S), Th2 (CD4.T/U), Th17 (CD4.X/Y) tip clusters",
    highlight_size = 1, highlight_alpha = 0.7,
    highlight_raster = TRUE, labelclusters = FALSE,
    which = "plot1"
))
dev.off()

# ============================================================
# S11b: Th-defining receptor transcripts on the CD4 MDE
# ============================================================
# The complete object stores counts only, so it is normalised here. Normalising
# after subsetting to CD4 is safe: the scale factor is the per-cell total over
# all genes, which does not depend on which cells are retained.
#
# Capped at the 95th percentile rather than rescaled per gene, so the colour
# scale is comparable across the six panels.
cd4 <- NormalizeData(cd4, assay = "RNA",
                     normalization.method = "LogNormalize", verbose = FALSE)

for (gene in TH_MARKERS) {
    panel_pdf(figure_dir, sprintf("S11b_CD4MDE_%s", gene), 5, 5)
    print(
        FeaturePlot(cd4, reduction = "mde_incremental", features = gene,
                    order = TRUE, cols = c("lightgrey", "blue"),
                    raster = TRUE, raster.dpi = c(512, 512),
                    max.cutoff = "q95") +
            ggtitle(gene)
    )
    dev.off()
}
