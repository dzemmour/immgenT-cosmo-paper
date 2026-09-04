# Figure 3. Reference MDE plot for each T-cell lineage.
#
# Panels produced:
#   3a  Lineage-specific MDE for each of the eight lineages, coloured and
#       labelled by cluster.
#   3b  CD44+CD62L- and CD62L+CD44- cells on each lineage-specific MDE, with
#       the protein gate they are drawn from.
#   3c  Dot plot of canonical T-cell transcripts across all clusters.
#
# --- internal ---
# Ported from cosmo_paper.Rmd sections "Level2 in lineage-specific MDE -
# Figure 3a" (line 1172), "CD62L CD44 plots Fig 3b" (line 1204) and "Figure 3c -
# Dot plot" (line 1259).
#
# The Rmd's "Fig 3b" section produced both this figure's panel b and Figure 2f;
# the two are distinguished by which embedding they use (mde_incremental,
# split by lineage, here; mde2_totalvi_20241006 for Figure 2f). It also gated
# double-positive and double-negative populations, which are not published.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   immgenT_seurat_ADT_GeneSubset.Rds   [primary input]

suppressPackageStartupMessages({
    library(dplyr)
    library(viridis)
})

source("code/R/setup.R")

figure_dir <- "Figure 3"
set.seed(1)

# The eight T-cell lineages, in the order the figure lays them out.
LINEAGES <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

# Point size per lineage: the lineage-specific MDEs hold between ~4,000 (DP) and
# ~200,000 (CD8) cells on the same physical panel size, so a single point size
# would leave the small lineages invisible and the large ones a solid block.
LINEAGE_PT_SIZE <- c(CD8 = 1, CD4 = 1, Treg = 2, gdT = 2,
                     CD8aa = 2, Tz = 2, DN = 3, DP = 3)

# ============================================================
# Load data
# ============================================================
so_orig <- readRDS(sprintf("%s/immgenT_seurat_ADT_GeneSubset.Rds", data_path))
so <- so_orig[, so_orig@meta.data %>%
                  filter(organ_simplified != "thymus") %>% rownames()]

# ============================================================
# 3a: Lineage-specific MDE coloured by cluster
# ============================================================
# Each lineage was embedded separately, so mde_incremental holds a different
# 2-D space for each; the object must therefore be split before plotting or
# cells from different lineages would be drawn on incompatible axes.
so_list <- SplitObject(so, split.by = "annotation_level1")

for (lin in LINEAGES) {
    panel_pdf(figure_dir, sprintf("3a_MDE_%s_level2", lin), 5, 5)
    print(
        DimPlot(
            so_list[[lin]], reduction = "mde_incremental",
            group.by = "annotation_level2", cols = mypal_level2,
            raster = TRUE, raster.dpi = c(1024, 1024),
            pt.size = LINEAGE_PT_SIZE[[lin]], alpha = 1,
            label = TRUE, label.size = 2, label.box = TRUE
        ) + NoGrid() + NoLegend()
    )
    dev.off()
}

rm(so_list)

# ============================================================
# 3b: CD62L/CD44 single-positive cells on each lineage MDE
# ============================================================
# Same gate as Figure 2f, applied here within each lineage. Density rather than
# individual points, because the question is where in each lineage the naive-
# and experienced-phenotype cells sit, not which individual cells they are.
so_cite <- so_orig[, so_orig@meta.data %>%
                       filter(organ_simplified != "thymus", cite_seq == TRUE) %>%
                       rownames()]
so_cite <- NormalizeData(so_cite, assay = "ADT",
                         normalization.method = "LogNormalize", verbose = FALSE)

cd62l_gate <- 4.5
cd44_gate <- 5

adt <- so_cite[["ADT"]]$data
so_cite$is_CD62L_sp <- adt["CD62L", ] > cd62l_gate & adt["CD44", ] < cd44_gate
so_cite$is_CD44_sp  <- adt["CD62L", ] < cd62l_gate & adt["CD44", ] > cd44_gate

so_cite$sub <- FALSE
so_cite$sub[sample(colnames(so_cite), size = 50000, replace = FALSE)] <- TRUE

panel_pdf(figure_dir, "3b_gate_CD62L_CD44", 5, 5)
print(
    MyFeatureScatter(
        so = so_cite[, so_cite$sub], assay = "ADT", slot = "data",
        feature1 = "CD62L", feature2 = "CD44", group.by = "sub",
        raster = TRUE, highlight_size = 1
    ) +
        geom_vline(xintercept = cd62l_gate, linetype = "dashed",
                   color = "brown") +
        geom_hline(yintercept = cd44_gate, linetype = "dashed",
                   color = "brown") +
        NoLegend()
)
dev.off()

panel_pdf(figure_dir, "3b_MDE_bylineage_CD44pos_CD62Lneg", 7, 7)
print(
    MyDimPlotHighlightDensity(
        seurat_object = so_cite, umap_to_plot = "mde_incremental",
        group.by = "is_CD44_sp", split.by = "annotation_level1",
        raster = TRUE, cols = plasma(10),
        highlight_size = 0.5, highlight_pixels = c(512, 512)
    ) + ggtitle("CD44+ CD62L-")
)
dev.off()

panel_pdf(figure_dir, "3b_MDE_bylineage_CD62Lpos_CD44neg", 7, 7)
print(
    MyDimPlotHighlightDensity(
        seurat_object = so_cite, umap_to_plot = "mde_incremental",
        group.by = "is_CD62L_sp", split.by = "annotation_level1",
        raster = TRUE, cols = plasma(10),
        highlight_size = 0.5, highlight_pixels = c(512, 512)
    ) + ggtitle("CD44- CD62L+")
)
dev.off()

# ============================================================
# 3c: Canonical transcripts across clusters
# ============================================================
# scale = FALSE: the dot colour is mean log1p CP10K expression, not a z-score
# across clusters. Scaling would make a gene expressed at a low level
# everywhere look as structured as a genuine lineage marker.
genes <- c("Cd3e", "Cd8a", "Cd8b1", "Cd4", "Trdc", "Zbtb16", "Foxp3", "Rorc",
           "Sell", "Cd44", "Itgae", "Ifng", "Il4", "Il17a", "Gzma", "Pdcd1",
           "Klrg1")

panel_pdf(figure_dir, "3c_DotPlot_canonical_genes", 25, 8)
print(
    Seurat::DotPlot(so_cite, features = rev(genes),
                    group.by = "annotation_level2", scale = FALSE) +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_viridis_c(option = "C")
)
dev.off()
