# Cluster composition helpers.
#
# Two related but distinct questions:
#   CompositionOfEachCluster()      -- which experiments contribute to a cluster
#                                      (Extended Data Figure 4a)
#   make_comp_dotdata_and_plot()    -- what proportion of a lineage each cluster
#                                      accounts for, per tissue or perturbation
#                                      (Figure 6b-d)
#   load_sample_composition()       -- the per-sample proportion table both
#                                      Figure 4d and Extended Data Figure 4b are
#                                      built from

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(rlang)
    library(tibble)
    library(reshape2)
    library(RColorBrewer)
})

#' Experiment composition of each cluster
#'
#' Counts cells per cluster x sample, then normalises twice: first within
#' sample, so that large samples do not dominate, and then within cluster, so
#' the result is the relative contribution of each sample to that cluster and
#' sums to 1. Without the first normalisation the plot would mostly report
#' which experiments sequenced the most cells.
#'
#' @param meta A metadata data frame.
#' @param cluster_col,sample_col Column names to cross-tabulate.
#' @return Long data frame with `p_in_sample`, `p_weighted` and a
#'   "cluster (n=XX)" display label ordered as the clusters appear.
CompositionOfEachCluster <- function(meta,
                                     cluster_col = "annotation_level2",
                                     sample_col  = "IGT") {
    clust  <- sym(cluster_col)
    sample <- sym(sample_col)

    comp <- meta %>%
        filter(!is.na(!!clust), !is.na(!!sample)) %>%
        count(!!clust, !!sample, name = "n") %>%
        group_by(!!sample) %>%
        mutate(p_in_sample = n / sum(n)) %>%
        ungroup() %>%
        group_by(!!clust) %>%
        mutate(p_weighted = p_in_sample / sum(p_in_sample)) %>%
        ungroup()

    cluster_sizes <- comp %>%
        group_by(!!clust) %>%
        summarise(n_total = sum(n), .groups = "drop") %>%
        mutate(
            clust_chr = as.character(!!clust),
            label = paste0(clust_chr, " (n=", n_total, ")")
        ) %>%
        select(-clust_chr)

    comp %>%
        left_join(cluster_sizes, by = rlang::as_name(clust)) %>%
        mutate(label = factor(label, levels = cluster_sizes$label))
}

#' Cluster proportion dot plot across a sample-level grouping
#'
#' Takes the per-sample proportion table, averages each cluster's
#' within-lineage proportion across the samples of each group (tissue,
#' perturbation, ...), and returns both the aggregated data and the dot plot.
#' Dot fill is the mean proportion and dot size the log10 number of cells, so
#' that a high proportion resting on very few cells is visibly distinguishable
#' from a well-sampled one. Rows are ordered by hierarchical clustering on the
#' proportion profiles, with "baseline"/"healthy" pinned to the top as the
#' reference row.
#'
#' @param data Per-sample proportion table (see `load_sample_composition`).
#' @param so Seurat object, used only for its lineage and cluster factor
#'   orderings so that panels agree with the rest of the paper.
#' @param condition_col Sample-level column to group by.
#' @return List with `data` (the aggregated table) and `plot`.
make_comp_dotdata_and_plot <- function(data, so,
                                       condition_col = "condition_broad") {
    cond_sym <- sym(condition_col)

    data_long <- data %>%
        select(!!cond_sym, IGTHT, matches("\\.(prop|ncells)$")) %>%
        pivot_longer(
            cols = matches("\\.(prop|ncells)$"),
            names_to = c("cluster", "metric"),
            names_pattern = "^(.*)\\.(prop|ncells)$",
            values_to = "value"
        ) %>%
        pivot_wider(names_from = metric, values_from = value)

    comp_mean <- data_long %>%
        group_by(!!cond_sym, cluster) %>%
        summarise(
            mean_prop    = mean(prop, na.rm = TRUE),
            total_ncells = sum(ncells, na.rm = TRUE),
            .groups = "drop"
        )

    cluster_map <- so@meta.data %>%
        as.data.frame() %>%
        distinct(annotation_level1, annotation_level2)

    comp_mean <- comp_mean %>%
        left_join(cluster_map, by = c("cluster" = "annotation_level2"))

    mat <- comp_mean %>%
        select(!!cond_sym, cluster, mean_prop) %>%
        pivot_wider(
            names_from  = cluster,
            values_from = mean_prop,
            values_fill = list(mean_prop = 0)
        ) %>%
        column_to_rownames(var = as_name(cond_sym)) %>%
        as.matrix()

    organ_order <- if (nrow(mat) > 1) {
        rownames(mat)[hclust(dist(mat))$order]
    } else {
        rownames(mat)
    }

    if ("baseline" %in% organ_order) {
        organ_order <- c("baseline", setdiff(organ_order, "baseline"))
    }
    if ("healthy" %in% organ_order) {
        organ_order <- c("healthy", setdiff(organ_order, "healthy"))
    }

    comp_mean <- comp_mean %>%
        mutate(!!condition_col := factor(!!cond_sym, levels = organ_order)) %>%
        mutate(
            annotation_level1 = factor(annotation_level1,
                                       levels = levels(so$annotation_level1)),
            cluster = factor(cluster, levels = levels(so$annotation_level2)),
            # +1 so that clusters with zero cells still receive a finite size.
            size_log10 = log10(total_ncells + 1)
        )

    p <- ggplot(comp_mean, aes(x = cluster, y = !!cond_sym)) +
        geom_point(
            aes(size = size_log10, fill = mean_prop),
            shape = 21, colour = "black", stroke = 0.2
        ) +
        scale_fill_gradientn(
            colours = RColorBrewer::brewer.pal(9, "Blues"),
            name = "Mean prop"
        ) +
        scale_size(range = c(1, 8), name = "log10(# cells)") +
        facet_grid(. ~ annotation_level1, scales = "free_x", space = "free_x") +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            panel.grid = element_blank(),
            strip.text = element_text(face = "bold")
        )

    list(data = comp_mean, plot = p)
}

#' Per-sample cluster proportion table joined to sample metadata
#'
#' Reads the cached per-sample proportion table (each cluster's share of its
#' own lineage within a sample) and joins the curated sample metadata. Organ
#' and perturbation are returned as factors in anatomical / reference-first
#' order rather than alphabetically, because that ordering is what every
#' composition panel uses.
#'
#' @param data_dir Directory holding the two input files.
#' @return A data frame, one row per sample (IGTHT).
load_sample_composition <- function(data_dir = "data") {
    prop_data <- readRDS(sprintf(
        "%s/annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds", data_dir
    )) %>%
        dplyr::select(IGTHT, ends_with("ncells"), ends_with("prop"))

    m <- read.table(
        sprintf("%s/Sample_metadata_david_20260107_v11.csv", data_dir),
        header = TRUE, sep = ","
    )

    # Renamed here rather than in the source file so that the file on disk stays
    # consistent with the sample database it is exported from.
    m$IGT_10xlane_id <- NULL
    m$organism_id <- NULL
    m$IGT <- m$immgent_IGT_10xlane_id
    m$immgent_IGT_10xlane_id <- NULL
    m$HT <- m$hashtag_number
    m$hashtag_number <- NULL
    m$sample_id <- NULL
    m$count_unique <- NULL
    m$sample_code.1 <- NULL

    mdata <- merge(prop_data, m, by = "IGTHT",
                   all.x = TRUE, all.y = FALSE, suffixes = c(".x", ""))

    mdata$organ_simplified <- factor(
        mdata$organ_simplified,
        levels = c("blood", "spleen", "LN", "SLO", "bone marrow", "thymus",
                   "lung", "liver", "peritoneal cavity",
                   "colon epi", "small intestine epi", "colon LP",
                   "small intestine LP", "skin",
                   "mammary gland", "uterus", "placenta", "prostate",
                   "submandibular gland", "kidney", "pancreas", "CNS",
                   "synovial fluid")
    )

    mdata$condition_detailed_simplified <-
        relevel(factor(mdata$condition_detailed_simplified), ref = "baseline")

    mdata %>% arrange(organ_simplified, condition_detailed_simplified)
}

#' Reshape the per-sample proportion table to one row per sample x cluster
#'
#' @param mdata Output of `load_sample_composition()`.
#' @param level2_levels Cluster factor levels, for consistent ordering.
#' @param keep_ncells Also return the cell count per sample x cluster, needed
#'   when the panel filters on a minimum cluster size.
#' @return Long data frame with `annotation_level1` / `annotation_level2`,
#'   thymocyte clusters removed.
sample_composition_long <- function(mdata, level2_levels,
                                    keep_ncells = FALSE) {
    lineage_of <- function(x) {
        out <- rep(NA_character_, length(x))
        for (lin in c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")) {
            out[grepl(paste0("^", lin, "\\."), x)] <- lin
        }
        out[grepl("thymocyte", x)] <- "thymocyte"
        factor(out, levels = c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz",
                               "DN", "DP", "thymocyte"))
    }

    if (keep_ncells) {
        df <- mdata %>%
            dplyr::select(IGTHT, organ_simplified, condition_detailed_simplified,
                          ends_with("prop"), ends_with(".ncells"))
        df_long <- df %>%
            pivot_longer(
                cols = matches("\\.(prop|ncells)$"),
                names_to = c("variable", ".value"),
                names_pattern = "^(.*)\\.(prop|ncells)$"
            )
        df_long$value <- df_long$prop
    } else {
        df <- mdata %>%
            dplyr::select(IGTHT, organ_simplified, condition_detailed_simplified,
                          ends_with("prop"))
        df_long <- reshape2::melt(df)
        df_long$variable <- gsub("\\.prop", "", df_long$variable)
    }

    df_long$annotation_level1 <- lineage_of(as.character(df_long$variable))
    df_long$IGTHT <- factor(df_long$IGTHT, levels = mdata$IGTHT)
    df_long$variable <- factor(as.character(df_long$variable),
                               levels = level2_levels)
    df_long$annotation_level2 <- df_long$variable

    df_long %>%
        filter(annotation_level1 != "thymocyte", !is.na(variable))
}
