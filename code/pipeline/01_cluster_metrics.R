# Per-experiment cluster robustness metrics.
#
# Produces the three cached files Figure 4b/4c and Extended Data Figure 4d
# read:
#
#   BioConservationSilhouette_IGT.Rds  silhouette width per de novo cluster,
#                                      before and after integration  (Fig 4b)
#   ClusterSeparation_RNA_IGT.Rds      separability per consensus cluster in
#                                      transcriptome space           (Fig 4c)
#   ClusterSeparation_ADT_IGT.Rds      the same in protein space     (ED 4d)
#
# Each metric needs an independent clustering or PCA of every experiment, so
# this is the slow step the figure scripts avoid by reading its output. Run it
# only to regenerate those files.
#
# The silhouette metric uses clusters found de novo within each experiment, so
# that it measures whether integration preserves the structure an experiment
# sees on its own. The separability metrics use the immgenT consensus clusters,
# so they measure whether the atlas partition is visible within each individual
# experiment. The two questions are different and the metrics are not
# comparable.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

source("code/R/separability.R")

data_path <- "data"

so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_complete.Rds", data_path))

# --- Figure 4b: silhouette before vs. after integration ----------------------
# The integrated space used here is the totalVI run fitted with the
# experiment-and-sample covariate removed, so that the comparison is not
# circular: the embedding being scored was not told which experiment each cell
# came from.
sil_summary_list <- list()
for (i in unique(so_orig$IGT)) {
    message(i)
    sil_summary_list[[i]] <- BioConservationSilhouette(
        so_orig[, so_orig$IGT == i],
        reduc_integrated = "totalvi_20241008_rmIGTsample"
    )
    sil_summary_list[[i]]$IGT <- i
}
saveRDS(sil_summary_list,
        sprintf("%s/BioConservationSilhouette_IGT.Rds", data_path))

# --- Figure 4c and Extended Data Figure 4d: cluster separability -------------
# Clusters with 10 or fewer cells in an experiment are skipped: a centroid
# estimated from a handful of cells is too noisy for the score to mean
# anything.
cluster_separability_by_experiment <- function(so, assay,
                                               normalization.method,
                                               nfeatures, npcs,
                                               features = NULL) {
    so$level2_int <- as.integer(so[["annotation_level2"]][, 1])
    DefaultAssay(so) <- assay
    dict_l2 <- so@meta.data %>% count(annotation_level2, level2_int)

    out <- list()
    for (i in unique(so$IGT)) {
        message(i)
        tmp <- if (is.null(features)) so[, so$IGT == i] else
            so[features, so$IGT == i]

        subcl <- tmp@meta.data %>%
            count(annotation_level2) %>%
            filter(n > 10) %>%
            pull(annotation_level2) %>%
            as.character()
        if (length(subcl) <= 1) next

        tmp <- tmp[, tmp$annotation_level2 %in% subcl]
        tmp <- tmp %>%
            NormalizeData(normalization.method = normalization.method,
                          verbose = FALSE) %>%
            FindVariableFeatures(selection.method = "vst",
                                 nfeatures = nfeatures, verbose = FALSE) %>%
            ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
            RunPCA(features = VariableFeatures(.), npcs = npcs,
                   verbose = FALSE)

        sep <- cluster_separability_centroid(
            emb = Embeddings(tmp, "pca")[, 1:10, drop = FALSE],
            clusters = as.integer(tmp[["level2_int"]][, 1]),
            b_mode = "mean"
        ) %>% as.data.frame()

        sep$annotation_level2 <-
            dict_l2$annotation_level2[match(sep$cluster, dict_l2$level2_int)]
        sep$IGT <- i
        out[[i]] <- sep
    }
    out
}

# Transcriptome space: the standard 2,000 variable genes and 50 components.
sep_rna <- cluster_separability_by_experiment(
    so_orig[, so_orig@meta.data %>%
                filter(organ_simplified != "thymus") %>% rownames()],
    assay = "RNA", normalization.method = "LogNormalize",
    nfeatures = 2000, npcs = 50
)
saveRDS(sep_rna, sprintf("%s/ClusterSeparation_RNA_IGT.Rds", data_path))

# Protein space: CLR normalisation, and only the antibodies that passed manual
# QC as good or intermediate -- a poorly performing antibody contributes noise
# that would depress every cluster's score equally.
so_cite <- so_orig[, so_orig@meta.data %>%
                       filter(organ_simplified != "thymus",
                              cite_seq == TRUE) %>% rownames()]
markers <- so_cite[["ADT"]]@meta.features %>%
    filter(IGT1only == FALSE,
           classification %in% c("good", "intermediate")) %>%
    rownames()

sep_adt <- cluster_separability_by_experiment(
    so_cite, assay = "ADT", normalization.method = "CLR",
    nfeatures = 50, npcs = 10, features = markers
)
saveRDS(sep_adt, sprintf("%s/ClusterSeparation_ADT_IGT.Rds", data_path))
