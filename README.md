# immgenT-cosmo-paper

Analysis code for *immgenT: A Comprehensive Reference of Convergent T-cell
States in the Mouse* ([bioRxiv 10.64898/2026.01.30.702892](https://www.biorxiv.org/content/10.64898/2026.01.30.702892v3)):
an atlas of 682,935 mouse T cells from 734 samples across 80 experiments,
profiled by single-cell RNA-seq, 128-plex CITE-seq and paired αβ TCR
sequencing, and organised into eight lineages and 107 reproducible clusters.

**Site: <https://dzemmour.github.io/immgenT-cosmo-paper/>** — one page per
published figure, showing for each panel the code that produced it, the panel
itself, and a caption.

## Layout

- `analysis/` — the [workflowr](https://workflowr.github.io/workflowr/) site
  source, one page per figure. Build with `workflowr::wflow_build()`; output
  goes to `docs/`.
- `script/` — one R script per published figure
  (`Rscript script/FigureN.R`), writing one PDF per panel into
  `figures/<figure>/`. See `script/README.md` for the figure → script map and
  the note on figure numbering.
- `code/R/` — shared plotting and data-loading helpers, sourced by the figure
  scripts. `code/pipeline/` — the slow upstream steps that produced the cached
  inputs, kept as provenance rather than as a rebuild. See `code/README.md`,
  which also carries the data provenance table.
- `figures/` — the panel set, one PDF per panel. This is what the scripts
  produce and what the site displays.
- `output/` — accompanying tables (specificity scores, QC summaries,
  per-study annotation rates).
- `docs/` — the built site.
- `data/` — **not tracked in git**; see below.

## Data

The inputs are ~3.9 GB of Seurat objects and cached intermediates, deposited on
Zenodo rather than in this repository. Download the archive and unpack it into
`data/` at the repository root, so that for example
`data/immgenT_seurat_ADT_GeneSubset.Rds` exists. Every file, what reads it and
where it comes from is listed in the data provenance table in
`code/README.md`.

Two conventions in that data are worth knowing before reading any script:

- **RNA in `immgenT_seurat_ADT_GeneSubset.Rds` is already log-normalised.** It
  had to be, because log-normalisation divides by each cell's total counts over
  all genes and the object is gene-subsetted. No figure script normalises its
  RNA assay again.
- **ADT is deliberately not normalised.** Each panel that plots protein
  normalises the exact set of cells it gates, so that a threshold means the
  same thing throughout a panel.

## Reproducing a figure

```
Rscript script/Figure4.R          # writes figures/Figure 4/*.pdf
Rscript script/render_site_assets.R   # converts panels to the site's PNGs
```

then `workflowr::wflow_build()` to rebuild the pages. The code on each page is
displayed, not re-executed at build time, so building the site does not require
the data.

## License

MIT (see `LICENSE`). The bundled JavaScript libraries under `docs/site_libs/`
are third-party and carry their own licences.

## Requirements

R 4.4 with Seurat 5, and [`ZemmourLib`](https://github.com/dzemmour/ZemmourLib)
— an internal package providing the immgenT colour palettes
(`immgent_colors`), the highlight and density plotting helpers
(`MyDimPlotHighlight`, `MyDimPlotHighlightDensity`, `MyFeatureScatter`) and the
pre-rendered lineage MDE backgrounds (`AddEmbeddingRasterBackground`) that
several panels draw over. The remaining dependencies are on CRAN and
Bioconductor: dplyr, tidyr, tibble, purrr, rlang, forcats, reshape2, ggplot2,
ggrepel, ggalluvial, scales, scattermore, pheatmap, ComplexHeatmap, circlize,
RColorBrewer, pals, viridis, cluster, RANN, magick.

`script/render_site_assets.R` additionally needs Ghostscript on `PATH`
(`brew install ghostscript`). Each page's exact session is recorded in the
session information block at the foot of the built page.
