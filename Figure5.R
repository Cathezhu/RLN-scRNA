suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ------------- Figure 5 --------------
# Plotting code for SPP1 repair macrophage panels.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

mac_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#cab2d6", "#fda111", "#80b1d3", "#8dd3c7")

mac_path <- file.path(datadir, "Figure5", "Macrophage.rds")
if (file.exists(mac_path)) {
  Macrophage <- readRDS(mac_path)
  Idents(Macrophage) <- "Mac"

  # Figure 5A-B: SPP1 expression across time and macrophage subtype
  if ("SPP1" %in% rownames(Macrophage)) {
    p <- VlnPlot(Macrophage, features = "SPP1", group.by = "orig.ident", pt.size = 0) +
      theme_classic() +
      labs(x = "", y = "SPP1 expression")
    save_panel("Figure5A_SPP1_by_time.pdf", p, width = 5, height = 4)

    p <- VlnPlot(Macrophage, features = "SPP1", group.by = "Mac", pt.size = 0, cols = mac_palette) +
      theme_classic() +
      labs(x = "", y = "SPP1 expression")
    save_panel("Figure5B_SPP1_by_macrophage_subtype.pdf", p, width = 5, height = 4)
  }
}

# Figure 5C: D3 CellChat heatmap from precomputed plot
plot_path <- file.path(datadir, "Figure5", "Figure5C_D3_CellChat_heatmap.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("Figure5C_D3_CellChat_heatmap.pdf", p, width = 6, height = 5)
}

# Figure 5D: SuperCell/metacell correlation panel
# Expected columns: gene_x, gene_y, cor, group or cell/metacell-level table for plotting.
cor_path <- file.path(datadir, "Figure5", "SuperCell_SPP1_correlation.csv")
if (file.exists(cor_path)) {
  cor_df <- read.csv(cor_path, check.names = FALSE)
  if (all(c("SPP1", "FABP4", "FN1") %in% colnames(cor_df))) {
    p <- ggplot(cor_df, aes(x = SPP1, y = FABP4, color = FN1)) +
      geom_point(size = 1.2, alpha = 0.8) +
      scale_color_viridis_c() +
      theme_classic() +
      labs(x = "SPP1", y = "FABP4", color = "FN1")
    save_panel("Figure5D_SuperCell_SPP1_FABP4_FN1_correlation.pdf", p, width = 5, height = 4)
  }
}

# Figure 5E-G: macrophage pseudotime/branch/gene dynamics from precomputed plot objects
for (panel in c("Figure5E_macrophage_pseudotime", "Figure5F_branch_heatmap", "Figure5G_pseudotime_gene_dynamics")) {
  plot_path <- file.path(datadir, "Figure5", paste0(panel, ".rds"))
  if (file.exists(plot_path)) {
    p <- readRDS(plot_path)
    save_panel(paste0(panel, ".pdf"), p, width = 6, height = 5)
  }
}

# Figure 5H: Mac_5 vs Mac_4 GSEA/Reactome enrichment barplot
gsea_path <- file.path(datadir, "Figure5", "Mac5_vs_Mac4_GSEA.csv")
if (file.exists(gsea_path)) {
  gsea <- read.csv(gsea_path, check.names = FALSE)
  p <- ggplot(gsea, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027") +
    theme_classic() +
    labs(x = "", y = "NES")
  save_panel("Figure5H_Mac5_vs_Mac4_GSEA.pdf", p, width = 7, height = 5)
}
