# Extended Data Figure 8. Cluster distributions across tissues.
#
# Panel produced:
#   S8  All-T MDE at baseline for the five sparsely sampled sites, coloured by
#       lineage.
#
# --- internal ---
# The same plot as Figure 6a, drawn for the five organs that figure leaves out.
# Both organ lists and the plotting code live in code/R/organ_mde.R; the Rmd
# produced them in a single loop over every baseline organ (section "Fig 6a -
# S6 - MDE plots organ", line 1663, whose "S6" label is stale).
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source("code/R/setup.R")
source("code/R/organ_mde.R")

figure_dir <- "Extended Data Figure 8"
set.seed(1)

so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so_baseline <- baseline_cells(so_orig)

for (organ in ORGANS_FIGURES8) {
    plot_organ_mde(so_baseline, organ, figure_dir, "S8_MDE_baseline")
}
