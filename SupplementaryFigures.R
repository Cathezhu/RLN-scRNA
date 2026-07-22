suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ------------- Supplementary Figures --------------
# Plotting code for supplementary single-cell panels.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

global_palette <- c(
  "#E64B35FF", "#f4c7b4", "#F39B7FFF",
  "#fae6b8", "#dfad80", "#B09C85FF",
  "#91D1C2FF", "#00A087FF", "#8491B4FF",
  "#ae93b8", "#4DBBD5FF", "#3C5488FF",
  "#7E6148FF", "#98c98b"
)

# Fig S1A-D: global clusters, QC metrics, marker UMAPs
obj_path <- file.path(datadir, "SupplementaryFigures", "new_celltype.rds")
if (file.exists(obj_path)) {
  obj <- readRDS(obj_path)
  pal_use <- rep(global_palette, length.out = length(levels(factor(obj$new_celltype))))

  p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = TRUE) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  save_panel("FigS1A_unannotated_clusters_umap.pdf", p, width = 5, height = 4)

  qc_features <- intersect(c("nFeature_RNA", "percent.mt"), colnames(obj@meta.data))
  if (length(qc_features) > 0) {
    p <- VlnPlot(obj, features = qc_features, group.by = "new_celltype", pt.size = 0, cols = pal_use, ncol = length(qc_features)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_panel("FigS1B-C_QC_metrics.pdf", p, width = 10, height = 4)
  }

  marker_features <- intersect(c("CD68", "CD3E", "PDGFRA", "COL15A1", "ACTA2", "MPZ"), rownames(obj))
  if (length(marker_features) > 0) {
    p <- FeaturePlot(obj, features = marker_features, reduction = "umap", order = TRUE, ncol = 3) +
      theme(legend.position = "right")
    save_panel("FigS1D_major_marker_featureplots.pdf", p, width = 9, height = 6)
  }
}

# Fig S2A: GO dotplot for upregulated genes over time
plot_path <- file.path(datadir, "SupplementaryFigures", "FigS2A_GO_dotplot.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("FigS2A_GO_dotplot.pdf", p, width = 8, height = 5)
}

# Fig S2B-C: ribosomal/TLR gene dynamics
gene_dyn_path <- file.path(datadir, "SupplementaryFigures", "FigS2_gene_dynamics.csv")
if (file.exists(gene_dyn_path)) {
  dyn <- read.csv(gene_dyn_path, check.names = FALSE)
  p <- ggplot(dyn, aes(x = time, y = expression, group = gene, color = gene)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~panel, scales = "free_y") +
    theme_classic() +
    labs(x = "", y = "Expression")
  save_panel("FigS2B-C_gene_dynamics.pdf", p, width = 8, height = 5)
}

# Fig S3: macrophage volcano/GSEA/PROGENy/GO/SCENIC panels
for (panel in c("FigS3A_mac_volcano", "FigS3B_IFN_GSEA", "FigS3C_PROGENy", "FigS3D_GO_dotplot", "FigS3E_SCENIC")) {
  plot_path <- file.path(datadir, "SupplementaryFigures", paste0(panel, ".rds"))
  if (file.exists(plot_path)) {
    p <- readRDS(plot_path)
    save_panel(paste0(panel, ".pdf"), p, width = 7, height = 5)
  }
}

# Fig S4: T-cell CellChat, markers, G2M, proportions, PRF1
for (panel in c("FigS4A_CellChat_circle", "FigS4B_T_marker_heatmap", "FigS4C_G2M", "FigS4D_T_fraction", "FigS4E_PRF1")) {
  plot_path <- file.path(datadir, "SupplementaryFigures", paste0(panel, ".rds"))
  if (file.exists(plot_path)) {
    p <- readRDS(plot_path)
    save_panel(paste0(panel, ".pdf"), p, width = 7, height = 5)
  }
}

# Fig S5: fibroblast monocle3, ECM/SPP1, rankNet, macrophage Slingshot
for (panel in c("FigS5A_fib_monocle3_umap", "FigS5B_fib_monocle3_trajectory", "FigS5C_ECM_score", "FigS5D_SPP1_overall", "FigS5E_rankNet", "FigS5F_mac_slingshot", "FigS5G_mac_slingshot", "FigS5H_mac_slingshot")) {
  plot_path <- file.path(datadir, "SupplementaryFigures", paste0(panel, ".rds"))
  if (file.exists(plot_path)) {
    p <- readRDS(plot_path)
    save_panel(paste0(panel, ".pdf"), p, width = 7, height = 5)
  }
}
