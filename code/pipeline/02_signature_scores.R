# Per-cell signature scores.
#
# Produces the two cached score files that Extended Data Figures 5 and 11c
# read:
#
#   signature_stress_scores.csv  14 stress-related signatures   (ED Fig 5)
#   treg_core_score.txt          the core Treg signature        (ED Fig 11c)
#
# Both use Seurat's AddModuleScore, which subtracts the mean expression of a
# set of control genes matched to the signature genes' expression bins. That
# control step is what makes the scores comparable between signatures of
# different size and expression level, and it is also what makes them
# expensive: scoring fourteen signatures against 683,000 cells is the reason
# these are cached rather than recomputed per figure.
#
# Scores are computed on the complete object, before the gene subsetting of
# 00_build_gene_subset.R, so that the control bins are drawn from the whole
# transcriptome and signature genes are not silently missing.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

data_path <- "data"

so <- readRDS(sprintf("%s/immgenT_seurat_ADT_complete.Rds", data_path))
so <- NormalizeData(so, assay = "RNA",
                    normalization.method = "LogNormalize")

# --- Extended Data Figure 5: stress signatures -------------------------------
sig_list <- readRDS(sprintf("%s/signatures_stress.Rds", data_path))
sig_list_use <- lapply(sig_list, function(x) intersect(x, rownames(so)))

# seed fixed so the control-gene sampling, and therefore the scores, are
# reproducible.
so <- AddModuleScore(
    so,
    features = sig_list_use,
    name = paste0(names(sig_list_use), "_"),
    assay = "RNA", slot = "data",
    ctrl = 100, nbin = 24, seed = 1
)

# AddModuleScore appends an index to each name it is given ("ier_1"); strip it
# so the columns are named after the signatures.
colnames(so@meta.data) <- sub("_[0-9]+$", "", colnames(so@meta.data))

write.csv(
    so@meta.data %>% dplyr::select(cellID, all_of(names(sig_list))),
    file = sprintf("%s/signature_stress_scores.csv", data_path),
    row.names = FALSE, quote = FALSE
)

# --- Extended Data Figure 11c: core Treg signature ---------------------------
# The eleven-gene core Treg signature of Zemmour et al.
treg_core_genes <- intersect(
    c("Foxp3", "Ikzf2", "Il2ra", "Ctla4", "Il2rb", "Capg", "Hopx",
      "Tnfrsf4", "Tnfrsf18", "Tnfrsf9", "Izumo1r"),
    rownames(so)
)

so <- AddModuleScore(so, features = list(treg_core_genes), name = "treg_core",
                     assay = "RNA", slot = "data")

write.table(so$treg_core1,
            file = sprintf("%s/treg_core_score.txt", data_path),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
