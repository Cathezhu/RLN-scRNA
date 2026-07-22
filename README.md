# RLN-scRNA

This repository contains plotting code for the manuscript figures. The code is organized as one R script per main figure, plus one script for supplementary figures.

## Files

- `Figure1.R`: global atlas UMAP, composition, marker feature plots, phenotype heatmap.
- `Figure2.R`: time-course gene dynamics, temporal heatmaps, GO/module score panels.
- `Figure3.R`: macrophage and T-cell plotting panels, SCENIC and NicheNet visualization.
- `Figure4.R`: fibroblast/myofibroblast plotting panels, trajectory and CellChat panels.
- `Figure5.R`: SPP1/macrophage repair panels, SuperCell correlation, pseudotime and enrichment panels.
- `SupplementaryFigures.R`: supplementary single-cell figure panels.

## Data policy

The scripts assume processed objects or intermediate plotting tables have already been generated. Raw 10X loading, scVI integration, re-clustering, CellChat inference, NicheNet inference, monocle/Slingshot model fitting, SCENIC, and other heavy analyses are not included here.

Place processed inputs under `./data` or set:

```r
Sys.setenv(RLN_FIGURE_DATA = "path/to/processed/data")
Sys.setenv(RLN_FIGURE_OUTPUT = "path/to/output/figures")
```

Each script writes panel-level outputs such as `Figure1E_umap.pdf` or `Figure5H_Mac5_vs_Mac4_GSEA.pdf`.
