# RLN-scRNA

This repository contains plotting code for the manuscript figures. The code is organized as one R script per main figure; related supplementary figure panels are placed at the end of the corresponding main-figure script.

## Files

- `Figure1.R`: Figure 1 and Fig S1 global atlas UMAP, composition, marker feature plots, phenotype heatmap and QC panels.
- `Figure2.R`: Figure 2 and Fig S2 temporal heatmaps, GO dotplot, ribosomal percentage violin plot, and module score panels.
- `Figure3.R`: Figure 3 plus Fig S3/S4 macrophage and T-cell plotting panels, SCENIC and NicheNet visualization.
- `Figure4.R`: Figure 4 plus Fig S5A-E fibroblast/myofibroblast plotting panels, trajectory and CellChat panels.
- `Figure5.R`: Figure 5 plus Fig S5F-H SPP1/macrophage repair panels, SuperCell correlation, pseudotime and enrichment panels.
- `CellChat.R`: standalone CellChat inference from a named time-point CellChat list, plus SPP1 and differential-interaction visualizations.
- `rankNet_ordered.R`: custom ordered CellChat rankNet helper used by Figure 4I.
- `supercell_2_Seurat.R`: SuperCell-to-Seurat conversion helper used by Figure 5D.

## Data policy

The figure scripts assume processed objects or intermediate plotting tables have already been generated. Raw 10X loading, scVI integration, re-clustering, NicheNet inference, SCENIC, and other heavy analyses are not included here. `CellChat.R` is provided separately for the CellChat inference workflow.

Place processed inputs under `./data` or set:

```r
Sys.setenv(RLN_FIGURE_DATA = "path/to/processed/data")
Sys.setenv(RLN_FIGURE_OUTPUT = "path/to/output/figures")
```

Each script writes panel-level outputs such as `Figure1E.pdf`, `FigS1A.pdf`, or `Figure5H.pdf`.
