# Extended Data Figure 1. Quality control metrics.
#
# Panels produced:
#   S1a  Per-experiment percentage of cells retained after QC, and the
#        percentage removed by each individual QC criterion.
#   S1b  Cells retained per experiment and per sample (hashtag).
#
# --- internal ---
# Ported from cosmo_paper.Rmd section "Figure S1 - QC" (line 869). That section
# also wrote the two QC summary tables; the sample-level one is published as
# Extended Data Table 3 and is written to output/FigureS1/ here.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md:
#   Sample_metadata_withQC_IGT1-96.csv  [primary input]
#   color_palette_igt.csv               [curated input]

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(reshape2)
})

source("code/R/setup.R")
source("code/R/utils.R")

figure_dir <- "Extended Data Figure 1"
out_dir <- "output/FigureS1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The 80 experiments that passed overall QC and contribute to the atlas.
# Listed explicitly rather than derived, because the QC metadata file also
# carries experiments that were run but excluded, and the figure reports the
# experiments the atlas is built from.
IGT_INCLUDED <- c(
    "IGT1", "IGT2", "IGT3", "IGT5", "IGT6", "IGT7", "IGT8", "IGT9",
    "IGT10", "IGT11", "IGT12", "IGT13", "IGT14", "IGT15", "IGT16", "IGT17",
    "IGT18", "IGT19", "IGT20", "IGT21", "IGT22", "IGT23", "IGT24", "IGT25",
    "IGT26", "IGT27", "IGT28", "IGT29", "IGT30", "IGT31", "IGT32", "IGT33",
    "IGT35", "IGT36", "IGT37", "IGT38", "IGT39", "IGT40", "IGT43", "IGT44",
    "IGT45", "IGT46", "IGT47", "IGT48", "IGT49", "IGT50", "IGT51", "IGT52",
    "IGT54", "IGT55", "IGT56", "IGT58", "IGT59", "IGT60", "IGT61", "IGT64",
    "IGT65", "IGT66", "IGT67", "IGT68", "IGT69", "IGT70", "IGT73", "IGT74",
    "IGT75", "IGT76", "IGT77", "IGT78", "IGT79", "IGT80", "IGT81", "IGT82",
    "IGT85", "IGT87", "IGT88", "IGT90", "IGT91", "IGT92", "IGT95", "IGT96"
)

# ============================================================
# Load data
# ============================================================
# fileEncoding strips the UTF-8 byte-order mark this export carries; without
# it the first column name comes through with the BOM attached and the
# experiment identifier column cannot be found.
sample_qc <- read.csv(sprintf("%s/Sample_metadata_withQC_IGT1-96.csv",
                              data_path),
                      fileEncoding = "UTF-8-BOM")
sample_qc$IGT <- sample_qc$immgent_IGT_10xlane_id
sample_qc <- sample_qc %>% filter(IGT %in% IGT_INCLUDED)
sample_qc$IGT <- ReorderIGT(sample_qc$IGT)

mypal_igt <- load_igt_palette()

# ============================================================
# S1a: QC pass and loss rates per experiment
# ============================================================
# Sample-level counts are summed to the experiment, then expressed as a
# percentage of the pre-QC count: experiments differ by an order of magnitude
# in size, so raw counts are not comparable between them.
igt_qc <- sample_qc %>%
    group_by(IGT) %>%
    summarize(across(c(ncells_preqc, ncells_outliers_nGenes,
                       ncells_outliers_deadcells,
                       ncells_outliers_lowCountADT,
                       ncells_outliers_autofluorescence,
                       ncells_nonTcells, ncells_postqc), sum),
              .groups = "drop") %>%
    mutate(
        percent_postqc = ncells_postqc / ncells_preqc * 100,
        percent_RNA_nGenes = ncells_outliers_nGenes / ncells_preqc * 100,
        percent_RNA_deadcells = ncells_outliers_deadcells / ncells_preqc * 100,
        percent_Protein_lowCount =
            ncells_outliers_lowCountADT / ncells_preqc * 100,
        percent_Protein_nonspecific_signal =
            ncells_outliers_autofluorescence / ncells_preqc * 100,
        percent_nonTcells = ncells_nonTcells / ncells_preqc * 100
    )

write.csv(igt_qc, file.path(out_dir, "S1a_QC_summary_by_IGT.csv"),
          row.names = FALSE)
write.csv(sample_qc, file.path(out_dir, "S1_QC_summary_by_sample.csv"),
          row.names = FALSE)

vars_keep <- c("percent_postqc", "percent_RNA_nGenes",
               "percent_RNA_deadcells", "percent_Protein_lowCount",
               "percent_Protein_nonspecific_signal", "percent_nonTcells")

# Only the experiments outside the 1.5 x IQR whiskers are labelled; labelling
# all 80 would bury the plot.
df <- melt(igt_qc, id.vars = "IGT") %>%
    filter(variable %in% vars_keep) %>%
    group_by(variable) %>%
    mutate(
        q1 = quantile(value, 0.25, na.rm = TRUE),
        q3 = quantile(value, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        is_outlier = value < (q1 - 1.5 * iqr) | value > (q3 + 1.5 * iqr)
    ) %>%
    ungroup()

p_s1a <- ggplot(df, aes(x = variable, y = value)) +
    geom_boxplot() +
    geom_jitter(aes(color = IGT), width = 0, size = 2) +
    geom_text_repel(
        data = df %>% filter(is_outlier),
        aes(label = IGT, color = IGT), max.overlaps = 30
    ) +
    scale_color_manual(values = mypal_igt) +
    labs(x = NULL, y = "% of pre-QC cells") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 75, vjust = 0.5, size = 10)) +
    ggtitle("QC RNA and Protein by IGT") +
    NoLegend() + NoGrid()

panel_pdf(figure_dir, "S1a_QC_RNA_ADT_by_IGT", 10, 8)
print(p_s1a)
dev.off()

# ============================================================
# S1b: Cells retained per experiment and per sample
# ============================================================
# One point per sample, coloured by hashtag, plus an "All" column pooling every
# sample into a single boxplot.
sample_qc <- sample_qc %>%
    mutate(percent_postqc = ncells_postqc / ncells_preqc * 100)

sample_qc_melt <- melt(sample_qc, id.vars = c("IGT", "hashtag_number"))

p_s1b <- ggplot(
    sample_qc_melt %>% filter(variable == "ncells_postqc"),
    aes(y = as.numeric(value))
) +
    geom_boxplot(aes(x = "All"), outlier.shape = NA, fill = "grey",
                 alpha = 0.5) +
    geom_jitter(aes(x = IGT, color = hashtag_number), width = 0, size = 2) +
    scale_color_manual(values = mypal) +
    scale_x_discrete(limits = c("All", levels(sample_qc$IGT))) +
    labs(x = NULL, y = "# cells after QC") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 75, vjust = 0.5, size = 10)) +
    ggtitle("Cells after QC by IGT and sample")

panel_pdf(figure_dir, "S1b_ncells_by_IGT_and_sample", 12, 8)
print(p_s1b)
dev.off()
