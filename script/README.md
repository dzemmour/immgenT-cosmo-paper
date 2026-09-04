# script/

One R script per published figure. Each script:

- reads only from `data/` (never writes there),
- writes its panels as one PDF per panel into `figures/<figure>/`, under the
  numbering of the **published** paper,
- writes any accompanying table into `output/<script>/`,
- sources shared plotting and data-loading helpers from `code/R/` rather than
  redefining them.

Run any script from the repository root:

```
Rscript script/Figure4.R
```

## Figure → script map

| Figure | Script | Panels produced | Panels not reproduced |
|---|---|---|---|
| 1 | `Figure1.R` | 1b, 1d, 1e, 1f | 1a, 1c (schematics) |
| 2 | `Figure2.R` | 2a, 2b, 2c, 2e, 2f, 2h | 2d, 2g (TCR analysis) |
| 3 | `Figure3.R` | 3a, 3b, 3c | — |
| 4 | `Figure4.R` | 4a–4j | — |
| 5 | `Figure5.R` | 5e | 5a–5d, 5f |
| 6 | `Figure6.R` | 6a, 6b, 6c, 6d | — |
| 7 | `Figure7.R` | 7a, 7b, 7c, 7d | 7e–7g |
| 8 | `Figure8.R` | 8a, 8b, 8c–8m | — |
| 9 | *(out of scope)* | — | all (gene programs; see below) |
| 10 | `Figure10.R` | 10a, 10b, 10c, 10f | 10d, 10e |
| 11 | *(not a plot)* | — | manuscript table |
| S1 | `FigureS1.R` | S1a, S1b | — |
| S2 | *(out of scope)* | — | all (CITE-seq QC; see below) |
| S3 | `FigureS3.R` | S3a–S3f | — |
| S4 | `FigureS4.R` | S4a, S4b, S4d | S4c |
| S5 | `FigureS5.R` | S5a–S5k | — |
| S6 | `FigureS6.R` | S6b–S6h | S6a |
| S7 | `FigureS7.R` | S7a–S7d | — |
| S8 | `FigureS8.R` | S8 | — |
| S9 | `FigureS9.R` | S9a, S9b | — |
| S10 | `FigureS10.R` | S10a–S10f | — |
| S11 | `FigureS11.R` | S11a, S11b, S11c | — |

`S*` here means Extended Data Figure, following the file-naming convention of
the [immgenT-GP-analysis](https://github.com/AgueroZZ/immgenT-GP-analysis)
repository; the site pages are titled "Extended Data Figure N".

`code/README.md` gives the reason each unreproduced panel is missing.

**Figure 9** (T-cell states are defined by combinations of gene programs) is
produced by the gene-program factorization analysis and lives in its own
repository: <https://github.com/AgueroZZ/immgenT-GP-analysis>.
**Extended Data Figure 2** (CITE-seq quality control) belongs to the companion
CITE-seq manuscript.

## A note on the figure numbering

The single analysis file this repository was ported from (`cosmo_paper.Rmd`,
kept alongside the repository but not tracked in it) labelled its sections on
an earlier draft's numbering. Three differences matter when comparing the two:

- effector molecules are **Figure 7** here and were labelled Figure 8 there;
- transcription factors are **Figure 8** here and were labelled Figure 7 there;
- the tissue-resident memory / CD8.Q analysis is **Figure 10** here and was
  labelled Figure 9 there.

Two of its supplementary labels were also stale: the organ MDE gallery's
supplementary counterpart is Extended Data Figure 8, not "S6", and its
Extended Data Figure 3 and Figure 2 panels were produced by the same chunks.
Each script's header records where its panels came from.

## Regenerating the site's images

The pages under `analysis/` display PNGs converted from the panel PDFs. After
re-running any figure script:

```
Rscript script/render_site_assets.R
```

This converts every `figures/**/*.pdf` to
`analysis/assets/<Figure>/<panel>.png` with Ghostscript, mirrors the result to
`docs/assets/` (which `workflowr::wflow_build()` does not refresh on its own),
and fails if a page references a PNG that was not produced. Converting rather
than re-plotting means a panel whose layout is not fully determined by the data
— ggrepel labels, jitter, the subsampling several panels use — cannot render
differently on the site than in the published PDF.
