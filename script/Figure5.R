# Figure 5. T-RBI integrates external datasets into the immgenT framework and
# supports near-saturation of T-cell states.
#
# Panels produced:
#   5e  Gallery of the external studies mapped onto the all-T MDE, one plot per
#       study, query cells coloured by lineage over the atlas in grey. The cell
#       count and percentage annotated at cluster level printed on each plot are
#       written to output/Figure5/5e_dataset_annotation_rates.csv.
#
# Panels 5a-5d (the Miller et al. worked example: author tSNE, T-RBI
# projection, immgenT annotation, confidence scores) and 5f (discovery score
# distribution across studies, excluding non-T cells) are not reproduced here;
# see analysis/Figure5.Rmd and code/README.md.
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure 5 T-RBI" (line 2077), chunks
# "MDE per dataset with immgent in background" and the all-dataset overview
# above it.
#
# The Rmd's discovery-score chunk in this section reads
# TRBI_discovery_scores_tbl_merged_DatasetPCA_withnonT.Rds, i.e. the score table
# that retains non-T cells as positive controls. That is Extended Data Figure
# 6g/6h, not Figure 5f, and it lives in script/FigureS6.R. The Figure 5f
# variant computed without non-T cells is not cached in data/ and no script
# here regenerates it.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   trbi_17studies_diet_merged.Rds      [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scattermore)
})

source("code/R/setup.R")
source("code/R/trbi_plots.R")

figure_dir <- "Figure 5"
out_dir <- "output/Figure5"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================
so_atlas <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so_merged_orig <- readRDS(sprintf("%s/trbi_17studies_diet_merged.Rds",
                                  data_path))

# Drop non-T cells and cells that failed T-cell filtering, leaving the 335,042
# query T cells the paper reports. The non-T cells are drawn separately by
# script/FigureS6.R.
so_merged <- so_merged_orig[, !so_merged_orig$level1_final %in%
                                c("nonT", "unclear")]
# "nonconv" was renamed Tz during revision; the palettes and every other figure
# use the new name.
so_merged$level1_final[so_merged$level1_final == "nonconv"] <- "Tz"

message("query T cells: ", ncol(so_merged))

# ============================================================
# 5e: Gallery of mapped external studies
# ============================================================
bkrg <- trbi_background(so_atlas, "mde2_totalvi_20241006")
dat_frg <- trbi_foreground_df(so_merged, "mde_incremental_allT")

# 299,764 of the 335,042 query T cells carry all-T MDE coordinates, spanning 16
# of the 17 dataset accessions in the object. GSE199563 was mapped against a
# lineage-specific reference only and so has no position in the all-T
# embedding; it is therefore absent from this gallery, which is why the figure
# shows sixteen plots.
message("query cells positioned in the all-T MDE: ", nrow(dat_frg))
plotted_datasets <- sort(unique(as.character(dat_frg$dataset)))

# All studies together first, then one plot per study.
panel_pdf(figure_dir, "5e_MDE_all_studies_level1", 5, 5)
print(
    bkrg +
        geom_scattermore(data = dat_frg,
                         aes(dim1, dim2, color = level1_final),
                         pointsize = 2, pixels = c(1024, 1024)) +
        scale_color_manual(values = mypal_level1) +
        trbi_theme() + NoLegend()
)
dev.off()

# Cluster-level annotation rate per study. A cell counts as annotated when its
# level2_final value names an actual cluster; the values that name a lineage
# instead ("CD4", "CD8", ...) are cells resolved only to lineage, and count as
# unannotated alongside "not classified".
unannotated_level2 <- c("not classified", "CD4", "CD8", "Treg", "gdT",
                        "nonconv", "Tz", "CD8aa", "DN", "DP", "thymocyte")

annotation_rates <- so_merged@meta.data %>%
    as.data.frame() %>%
    filter(dataset %in% plotted_datasets) %>%
    group_by(dataset) %>%
    summarise(
        n_cells = n(),
        pct_annotated_level2 =
            100 * mean(!level2_final %in% unannotated_level2),
        .groups = "drop"
    )

write.csv(annotation_rates,
          file.path(out_dir, "5e_dataset_annotation_rates.csv"),
          row.names = FALSE)

for (ds in plotted_datasets) {
    rate <- annotation_rates %>% filter(dataset == ds)
    panel_pdf(figure_dir, sprintf("5e_MDE_%s", ds), 5, 5)
    print(
        bkrg +
            geom_scattermore(data = dat_frg %>% filter(dataset == ds),
                             aes(dim1, dim2, color = level1_final),
                             pointsize = 2, pixels = c(1024, 1024)) +
            scale_color_manual(values = mypal_level1) +
            ggtitle(sprintf("%s\n%s cells, %.0f%% annotated", ds,
                            format(rate$n_cells, big.mark = ","),
                            rate$pct_annotated_level2)) +
            trbi_theme() + NoLegend()
    )
    dev.off()
}
