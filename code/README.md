# code/

Reusable R code backing `script/`. Nothing here writes to `figures/` — see
`../script/README.md` for how the figure scripts use it.

## R/ — sourced helpers

| File | Used by | Contents |
|---|---|---|
| `setup.R` | every script | `data_path`, the colour palettes (`mypal_level1`, `mypal_level2`, `mypal_organ`, `mypal_level2group`, `mypal`, `ColorRamp`), `load_igt_palette()`, the atlas loader `load_immgent()`, the panel device `panel_pdf()`, and `highlight_mde()`. |
| `utils.R` | Figure4, FigureS1 | `ReorderIGT()`: orders experiment identifiers numerically rather than as strings, so IGT10 does not sort before IGT5. |
| `organ_mde.R` | Figure6, FigureS8 | `ORGANS_FIGURE6` / `ORGANS_FIGURES8` (the two organ lists and their panel order), `baseline_cells()`, `plot_organ_mde()`. Shared so the two figures cannot disagree about which organ belongs to which. |
| `composition.R` | Figure4, Figure6, FigureS4 | `CompositionOfEachCluster()`, `make_comp_dotdata_and_plot()`, `load_sample_composition()`, `sample_composition_long()`. |
| `separability.R` | `pipeline/01` | `BioConservationSilhouette()`, `cluster_separability_centroid()`. Not called by any figure script — the figures read the cached metrics — but retained so those files have a documented producer. |
| `gene_auc.R` | Figure7, Figure8, FigureS9, FigureS10 | `gene_specificity_auc()` and the three plots built on it: `plot_gene_auc_curves()`, `plot_gene_barplot()`, `plot_gene_heatmap()` (with `heatmap_annotation()`). |
| `trbi_plots.R` | Figure5, FigureS6 | `trbi_background()`, `trbi_foreground_df()`, `trbi_theme()`, and the fixed MDE axis limits every T-RBI panel shares. |
| `discovery.R` | *(none — provenance only)* | `knn_novelty_scores()`, `summarize_discovery()`: the discovery score of Figure 5f and Extended Data Figure 6g/6h. The score table is cached and no script here regenerates it (see the gap noted below). |
| `alluvial.R` | Figure10 | `PlotAlluvialHighlight()`. |
| `doc_code.R` | `analysis/*.Rmd` | `code_for()`: selects a block of a figure script by *named anchor* so the site's code listings follow the code rather than going stale against hard-coded line numbers. Taken unchanged from the [immgenT-GP-analysis](https://github.com/AgueroZZ/immgenT-GP-analysis) repository. |

### `# --- internal ---` markers

`doc_code.R` strips anything between `# --- internal ---` and
`# --- end internal ---` from the code it puts on the site. The published pages
are read by people outside the project, so provenance notes that only make
sense internally — which section of the original analysis file a panel came
from, what an earlier draft numbered it, which line held a typo — live between
these markers: kept in the source for us, kept off the page.

## pipeline/ — data preparation, upstream of the figure scripts

Numbered in dependency order. **These are provenance, not a one-command
rebuild.** All three need the complete 2.7 GB atlas object and are slow enough
that their outputs are distributed as cached files; the figure scripts read
those files and never run these steps.

1. `00_build_gene_subset.R` — derives `immgenT_seurat_ADT_GeneSubset.Rds`, the
   object nearly every figure reads, from the complete object. Note the
   ordering constraint documented in its header: RNA normalisation must happen
   before gene subsetting, which is why the subsetted object ships
   already-normalised and no figure script normalises RNA again.
2. `01_cluster_metrics.R` — the three per-experiment cluster robustness metrics
   behind Figure 4b, Figure 4c and Extended Data Figure 4d.
3. `02_signature_scores.R` — the per-cell module scores behind Extended Data
   Figure 5 and Extended Data Figure 11c.

## Data provenance

Every `data/` file read anywhere in `script/`, and which step (if any) in this
repository produces it. `data/` is not tracked in git — it is ~3.9 GB and is
distributed via Zenodo as `immgenT-Cosmo.zip`
(<https://doi.org/10.5281/zenodo.21839963>), to be unpacked at the repository
root. Figure scripts also flag their own inputs inline via a `Required inputs`
header comment pointing back at this table.

| `data/` file | Produced by | Notes |
|---|---|---|
| `immgenT_seurat_ADT_complete.Rds` | *(primary input)* | The processed atlas: 682,935 cells × 55,494 genes, RNA counts plus 128-plex ADT, with the totalVI latent spaces and both MDE embeddings. The starting point; not produced by anything here. Read directly only by `FigureS11.R` and by `pipeline/`. |
| `immgenT_seurat_ADT_GeneSubset.Rds` | `00_build_gene_subset.R` | The same object cut to 2,058 genes and 128 proteins, 589 MB. **RNA is already log-normalised**; ADT is not. Read by every figure script except `FigureS1.R`. |
| `color_palette_igt.csv` | *(curated input)* | Fixed colour per experiment, so an experiment keeps its colour between panels that show different subsets of experiments. |
| `IGT5_seurat.Rds`, `IGT20_seurat.Rds`, `IGT21_seurat.Rds`, `IGT27_seurat.Rds`, `IGT40_seurat.Rds` | *(primary inputs)* | Per-experiment objects as delivered by the per-experiment pipeline, with each experiment's own UMAPs and de novo clusters. Used by the before/after-integration panels (Figures 1, 2h, 4f, 4i). Seurat 4.1.3 objects carrying `qc_stats_*` slots that hold plain data frames; see `script/Figure1.R` for why metadata is written through `@meta.data` on these. |
| `Sample_metadata_withQC_IGT1-96.csv` | *(primary input)* | Per-sample QC counts. Carries a UTF-8 byte-order mark, which `FigureS1.R` strips explicitly. Published as Extended Data Table 3. |
| `Sample_metadata_david_20260107_v11.csv` | *(curated input)* | Curated sample-level metadata, joined to the composition table by `load_sample_composition()`. Published as Extended Data Table 1. |
| `annotation_level2_PropPerSamplePerLevel1_withMetadata.Rds` | *(primary input)* | Per-sample proportion of each cluster **within its lineage**, plus cell counts. The within-lineage denominator is what makes tissues of very different lineage composition comparable. Feeds Figure 4d, Figure 6b-d and Extended Data Figure 4b. |
| `BioConservationSilhouette_IGT.Rds` | `01_cluster_metrics.R` | Silhouette width per de novo cluster per experiment, pre- and post-integration (Figure 4b). |
| `ClusterSeparation_RNA_IGT.Rds` | `01_cluster_metrics.R` | Centroid separability per consensus cluster per experiment, transcriptome space (Figure 4c). |
| `ClusterSeparation_ADT_IGT.Rds` | `01_cluster_metrics.R` | The same in protein space, on manually QC-passed antibodies only (Extended Data Figure 4d). |
| `TF_list.txt` | *(curated input)* | Transcription factor list for Figure 8 and Extended Data Figure 10. |
| `effmol_curated.txt` | *(curated input)* | 59 curated effector genes; 54 pass the expression filter and are the rows of Figure 7a. |
| `effmol_GO0005615_curated.txt` | *(curated input)* | 1,821 genes annotated to GO:0005615 (extracellular space); 305 pass the filter and the 100 most cluster-specific are the rows of Extended Data Figure 9a. This is the only difference between that heatmap and Figure 7a. |
| `igt1_96_pseudobulk_byclusterannotationlevel1_log1pCP10K.Rds` | *(primary input)* | Pseudobulk log1p CP10K per lineage, 8 columns (Figure 8a, Extended Data Figure 10a). |
| `igt1_96_pseudobulk_byclusterannotationlevel2_log1pCP10K.Rds` | *(primary input)* | Pseudobulk log1p CP10K per cluster, 91 columns (Figures 7, 8, Extended Data Figures 9, 10). |
| `signatures_stress.Rds` | *(curated input)* | 14 stress-related gene signatures, drawn from the published sources named in the Extended Data Figure 5 legend. |
| `signature_stress_scores.csv` | `02_signature_scores.R` | Per-cell scores for those 14 signatures (Extended Data Figure 5). |
| `treg_core_score.txt` | `02_signature_scores.R` | Per-cell score for the 11-gene core Treg signature (Extended Data Figure 11c). |
| `trbi_17studies_diet_merged.Rds` | *(primary input)* | The external studies after T-RBI mapping: 338,243 cells, of which 335,042 are T cells, carrying level-1/level-2 assignments, confidence scores and anchored MDE coordinates. 16 of the 17 dataset accessions have all-T MDE coordinates (299,764 cells); GSE199563 was mapped against a lineage-specific reference only, which is why the Figure 5e gallery shows sixteen plots. |
| `TRBI_discovery_scores_tbl_merged_DatasetPCA_withnonT.Rds` | **gap** | Per-cell discovery score in each study's own PCA space, retaining non-T cells as positive controls (Extended Data Figure 6g/6h). `R/discovery.R` holds the scoring functions, but the driver that ran them across studies is not preserved; the original analysis file read this table without producing it. |
| `CD4Ablation_trbi_seurat_objects.Rds` | *(primary input)* | One reduced Seurat object per ablated CD4 cluster (18 in all), holding the assignment and confidence columns before and after fine-tuning plus the ablated-reference MDE. |
| `CD4NoAblation_trbi_seurat.Rds` | *(primary input)* | The same 20,277 CD4 query cells mapped against the complete reference, as the no-ablation comparison for Extended Data Figure 7c. |
| `CD4NoAblation_trbi_predictions.csv` | *(unused)* | Per-cell predictions from the no-ablation run. Superseded by the columns carried on the object above; no script reads it. |
| `QC_RNA_ADT_sample_summary_table.tsv`, `QC_ncells_perIGTandsample.pdf` | *(outputs of an earlier run)* | Earlier copies of what `FigureS1.R` now writes into `output/FigureS1/`. Kept for comparison; nothing reads them. |

### Panels not reproduced, and why

| Panel | Missing input |
|---|---|
| Figure 1a, 1c; Figure 9a | Schematics drawn outside R. |
| Figure 2d, 2g | Clonotype-level TCR tables, from the TCR analysis. |
| Figure 5a-5d | The Miller et al. study's own tSNE embedding and author cluster labels. |
| Figure 5f | The discovery-score table computed with non-T cells excluded; only the with-non-T variant is cached (see the gap above). |
| Figure 7e-7g | A pseudobulk matrix computed per sample; only per-lineage and per-cluster matrices are cached. |
| Figure 9 (all) | The gene-program factorization, documented in the [immgenT-GP-analysis](https://aguerozz.github.io/immgenT-GP-analysis/) repository. |
| Figure 10d, 10e | The Milner et al. Trm gene signature and that study's dataset. |
| Figure 11 | A table typeset for the manuscript, not a plot. |
| Extended Data Figure 2 (all) | The CITE-seq per-antibody evaluation, documented in the companion CITE-seq manuscript; results published as Extended Data Table 4. |
| Extended Data Figure 4c | The per-cluster cross-experiment cosine similarity matrix. |
| Extended Data Figure 6a | The Miller et al. study's own UMAP embedding. |
| Extended Data Figure 7d | The per-cell discovery score under ablation. The cached ablation objects were reduced to the assignment and confidence columns and do not retain it. |

None of these gaps blocks reproducing the panels that are listed as produced;
they only matter for the specific panels named.
