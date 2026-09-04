# Figure 1. Scope of the immgenT Open Source Project.
#
# Panels produced:
#   1b  Cell counts by organ (top bars) and by broad disease category (left
#       bars). The condition table between them is typeset in the figure, not
#       plotted.
#   1d  IGT5 transcriptome UMAP coloured by transcriptomic cluster.
#   1e  IGT5 surface-protein UMAP, coloured by the clusters of (d).
#   1f  IGT5 flow-cytometry-like protein plots: TCRgd vs TCRb, CD4 vs CD8b.
#
# Panels 1a and 1c are schematics drawn outside R and are not produced here.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Figure 1 - Counting cells in immgenT"
# (line 729) and "Figure 1c-f (ok)" (line 992). The Rmd's counting section also
# produced a per-organ barplot at full organ resolution, a by-sex barplot and a
# by-perturbation barplot; none of the three is a published panel and they are
# dropped here. Its IGT5 section additionally plotted the hashtag (sample) UMAP,
# which is likewise not a published panel.
#
# The Rmd recoded condition_broad "autoimmunity" -> "autoimmune", which is a
# no-op against the current object: condition_broad is already curated. Its
# second recode was spelled "mutiinfection" and so never fired, but the merge it
# intended is real -- the published Virus bar is virus + multiinfection -- so it
# is reproduced here under the correct spelling.
#
# cosmo_paper.Rmd line 994 held a stray `sprintf(, wd)` that is a syntax error;
# it did nothing and is not reproduced.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]
#   IGT5_seurat.Rds                     [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(forcats)
    library(scales)
})

source("code/R/setup.R")

figure_dir <- "Figure 1"

# ============================================================
# Load data
# ============================================================
so <- load_immgent()

# ============================================================
# 1b: Cell counts by organ and by broad disease category
# ============================================================
# Organs are collapsed to the seven groups the figure shows. Spleen, lymph
# node, other secondary lymphoid organs and blood are pooled as the circulating
# / lymphoid compartment; colon and small intestine as gut; and the remaining
# sparsely sampled sites as "other".
df <- so@meta.data
df$organ_grouped <- as.character(df$organ_simplified)
df$organ_grouped[df$organ_simplified %in%
                     c("spleen", "LN", "SLO", "blood")] <- "spleen LN blood"
df$organ_grouped[grepl("colon|small intestine", df$organ_simplified)] <- "gut"
df$organ_grouped[grepl(paste0("synovial|prostate|uterus|pancreas|peritoneal|",
                             "bone|kidney|liver|mammary|placenta|thymus"),
                       df$organ_simplified)] <- "other"

df_organ <- df %>%
    filter(!is.na(organ_grouped)) %>%
    count(organ_grouped, name = "n")
df_organ$organ_grouped <- factor(
    df_organ$organ_grouped,
    levels = c("spleen LN blood", "gut", "lung", "skin",
               "submandibular gland", "CNS", "other")
)

p_1b_organ <- ggplot(df_organ, aes(x = organ_grouped, y = n,
                                   fill = organ_grouped)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = comma(n)), vjust = -0.25, size = 3) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0.02, 0.10))) +
    labs(title = "# cells by organ", x = NULL, y = "cells") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
    ) +
    NoGrid()

panel_pdf(figure_dir, "1b_cellcount_organ", 10, 7)
print(p_1b_organ)
dev.off()

# Co-infection experiments are reported under "virus": every experiment in the
# multiinfection category pairs a viral challenge with a second agent, and the
# figure shows ten disease categories, not eleven. Merging here reproduces the
# published Virus bar exactly (116,922 + 11,676 = 128,598 cells).
df$condition_broad[df$condition_broad == "multiinfection"] <- "virus"

df_cond <- df %>%
    filter(!is.na(condition_broad)) %>%
    count(condition_broad, name = "n")
# Reversed because coord_flip() draws the first factor level at the bottom;
# this puts healthy at the top of the panel.
df_cond$condition_broad <- factor(
    df_cond$condition_broad,
    levels = rev(c("healthy", "allergy", "autoimmune", "bacteria", "virus",
                   "parasite", "fungus", "cancer", "allotransplant",
                   "immunization"))
)

# Bars run leftwards (-n) after the flip so the category labels sit on the
# right-hand side, matching the figure layout.
p_1b_cond <- ggplot(df_cond, aes(x = condition_broad, y = -n,
                                 fill = condition_broad)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = comma(n)), hjust = 0, size = 3) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0.02, 0.10))) +
    labs(title = "# cells by disease category", x = NULL, y = "cells") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
    ) +
    NoGrid()

panel_pdf(figure_dir, "1b_cellcount_conditionbroad", 5, 5)
print(p_1b_cond)
dev.off()

# ============================================================
# 1d-1f: Representative experiment IGT5 (skin, baseline)
# ============================================================
# IGT5 is shown as delivered by the per-experiment pipeline: its own
# transcriptome and protein UMAPs and its own de novo clusters, before any
# cross-experiment integration. That is the starting point the rest of the
# paper integrates away from.
igt <- readRDS(sprintf("%s/IGT5_seurat.Rds", data_path))

# Written through @meta.data rather than `igt$all_cells <-`. The
# per-experiment objects are Seurat 4.1.3 and carry several `qc_stats_*` slots
# that hold plain data frames rather than Assay objects; under Seurat 5 the
# `[[<-` method revalidates the whole object and rejects them, while a direct
# @meta.data assignment does not. The same applies to every per-experiment
# object used in Figures 2 and 4.
igt@meta.data$all_cells <- TRUE

panel_pdf(figure_dir, "1d_IGT5_UMAP_RNA_clusters", 5, 5)
print(DimPlot(igt, reduction = "umap_rna", group.by = "RNA_clusters",
              cols = glasbey()))
dev.off()

panel_pdf(figure_dir, "1e_IGT5_UMAP_ADT_RNAclusters", 5, 5)
print(DimPlot(igt, reduction = "umap_adt", group.by = "RNA_clusters",
              cols = glasbey()))
dev.off()

# --- doc:1f ---
# group.by is a constant TRUE column: MyFeatureScatter highlights the grouped
# cells over a density background, and here every cell is to be highlighted.
panel_pdf(figure_dir, "1f_IGT5_TCRB_TCRGD", 5, 5)
print(MyFeatureScatter(so = igt, assay = "ADT", slot = "data",
                       feature1 = "TCRB", feature2 = "TCRGD",
                       group.by = "all_cells", raster = FALSE,
                       highlight_size = 2))
dev.off()

panel_pdf(figure_dir, "1f_IGT5_CD4_CD8B", 5, 5)
print(MyFeatureScatter(so = igt, assay = "ADT", slot = "data",
                       feature1 = "CD4", feature2 = "CD8B",
                       group.by = "all_cells", raster = TRUE,
                       highlight_size = 2))
dev.off()
