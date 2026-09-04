# Build the gene- and protein-subsetted atlas object every figure script reads.
#
#   immgenT_seurat_ADT_complete.Rds  ->  immgenT_seurat_ADT_GeneSubset.Rds
#   2.7 GB, 55,494 genes                589 MB, 2,058 genes
#
# The complete object is too large to load repeatedly, and no figure in this
# paper needs more than a few thousand genes. This step keeps the genes the
# analysis actually uses and drops the rest.
#
# CRITICAL: RNA normalisation happens *before* subsetting. Log-normalisation
# divides each cell by its total counts over all genes, so normalising after
# the subset would rescale every value by a different, gene-set-dependent
# factor. The subsetted object therefore ships with an already-normalised
# `data` layer, and no figure script normalises its RNA assay again.
#
# The ADT assay is subsetted but deliberately left un-normalised: each panel
# that plots protein normalises the exact set of cells it gates, so that
# thresholds mean the same thing within a panel.
#
# Not run as part of any figure. Requires the complete object.

suppressPackageStartupMessages({
    library(Seurat)
})

data_path <- "data"

so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_complete.Rds", data_path))
so_orig <- NormalizeData(so_orig, assay = "RNA",
                         normalization.method = "LogNormalize")

# Genes excluded and why:
#   TCR V/D/J/C segments  -- clonotype-driven, would make the embedding track
#                            clonal identity rather than cell state
#   Gm*/*Rik/*-ps         -- unannotated loci and pseudogenes
#   ribosomal, mitochondrial -- dominated by library size and cell viability
#   immunoglobulin        -- contaminating B-cell transcripts
#   never detected        -- carry no information
message("Removing TCR, mitochondrial, ribosomal, Gm/Rik, Ig and undetected genes")
tcr_genes <- grepl(paste0("Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|",
                          "Trgc|Trdv|Trdj|Trdc"), rownames(so_orig))
gm_rik_genes <- grepl("Gm|Rik$|\\-ps$", rownames(so_orig))
ribo_genes <- grepl("Rpl|Rps|Mrpl|Mrps|Rsl", rownames(so_orig))
mt_genes <- grepl("^mt-", rownames(so_orig))
ig_genes <- grepl("^Igh|^Igk|^Igl|^Igha|^Ighm|^Ighg|^Ighd|^Ighe",
                  rownames(so_orig))
genes_not_expressed <- rowSums(so_orig[["RNA"]]$counts) == 0

genes_to_keep <- !(tcr_genes | gm_rik_genes | ribo_genes | mt_genes |
                       genes_not_expressed | ig_genes)
cat("Number of genes to keep:", sum(genes_to_keep), "\n")

# IGT1 was stained with a smaller antibody panel; proteins unique to it cannot
# be compared across experiments.
prot_to_keep <- rownames(so_orig[["ADT"]]@meta.features)[
    so_orig[["ADT"]]@meta.features$IGT1only == FALSE
]

# Genes retained on top of the filter: the curated transcription factor and
# effector lists (Figures 7-8), the markers named in Figures 1-3 and 10, and
# the 1,000 most variable genes, so the object supports exploratory work as
# well as the published panels.
tf <- read.table(sprintf("%s/TF_list.txt", data_path),
                 header = FALSE, stringsAsFactors = FALSE)[, 1]
effmol <- read.table(sprintf("%s/effmol_curated.txt", data_path),
                     header = TRUE)[, 1]
genes_sub <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1", "Trdc", "Zbtb16", "Foxp3",
               "Rorc", "Sell", "Cd44", "Itgae", "Ifng", "Il4", "Il17a",
               "Gzma", "Pdcd1", "Klrg1", "Mki67", "Cdkn2a", "Cdkn2b",
               "Ikzf2", "Il2ra", "Ctla4", "Il2rb", "Capg", "Hopx",
               "Tnfrsf4", "Tnfrsf18", "Tnfrsf9", "Izumo1r")

so_orig <- FindVariableFeatures(so_orig, selection.method = "vst",
                                nfeatures = 1000)
hvg_1000 <- VariableFeatures(so_orig)
hvg_1000 <- hvg_1000[hvg_1000 %in% rownames(so_orig)[genes_to_keep]]

features_keep <- unique(c(tf, effmol, genes_sub, hvg_1000))

so_diet <- DietSeurat(
    so_orig,
    features  = c(features_keep, prot_to_keep),
    assays    = c("RNA", "ADT"),
    layers    = c("counts", "data"),
    dimreducs = c("mde2_totalvi_20241006", "mde_incremental"),
    graphs    = NULL
)

saveRDS(so_diet, sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
