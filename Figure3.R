suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(clusterProfiler)
})

# ------------- Figure 3 --------------
# Plotting code for macrophage and T-cell panels.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

mac_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#cab2d6", "#fda111", "#80b1d3", "#8dd3c7")
t_palette <- c("#E47BC0", "#8C87B3", "#E3C4DA", "#89C4D2", "#91D1C2FF", "#89E1F7", "#DFDDF0", "#F27F96")

mac_path <- file.path(datadir, "Figure3", "Macrophage.rds")
if (file.exists(mac_path)) {
  Macrophage <- readRDS(mac_path)
  Idents(Macrophage) <- "Mac"

  # Figure 3A: macrophage UMAP
  p <- DimPlot(Macrophage, reduction = "umap", group.by = "Mac", cols = mac_palette[seq_len(nlevels(Macrophage))], label = TRUE, repel = TRUE, raster = TRUE) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  save_panel("Figure3A_macrophage_umap.pdf", p, width = 5, height = 4)

  # Figure 3B: macrophage marker dotplot
  genes <- list(
    Mac_1 = c("CCL14", "C1QA", "C1QB", "SELENOP", "LYVE1"),
    Mac_2 = c("CCR2", "FCN1", "MEFV", "SELL"),
    Mac_3 = c("MX2", "MAT2A", "SSH2", "DHX34", "ND2"),
    Mac_4 = c("FABP4", "LPL", "FABP5", "MMP19", "CD36"),
    Mac_5 = c("EGLN3", "FLT1", "SERPINE1", "PCOLCE2", "ARG1", "FN1"),
    Cycling = c("CCNA2", "PBK", "KIF2C", "TOP2A", "MKI67")
  )
  p <- DotPlot(Macrophage, features = genes, group.by = "Mac") +
    coord_flip() +
    scale_color_gradientn(colours = c("#1B4F72", "#F5F5F5", "#922B21"), limits = c(-2.5, 2.5), oob = scales::squish) +
    theme_bw() +
    labs(x = "", y = "")
  save_panel("Figure3B_macrophage_marker_dotplot.pdf", p, width = 6, height = 7)

  # Figure 3C: macrophage composition across time
  meta <- Macrophage@meta.data %>%
    count(orig.ident, Mac, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()
  p <- ggplot(meta, aes(x = orig.ident, y = Mac)) +
    geom_point(aes(size = percent, color = percent)) +
    scale_color_viridis_c(option = "magma") +
    theme_bw() +
    labs(x = "", y = "", color = "Percent", size = "Percent")
  save_panel("Figure3C_macrophage_composition.pdf", p, width = 5, height = 4)
}

# Figure 3D/G: module-score boxplots for macrophage functions
# Expected columns: score_name, Mac, score
score_path <- file.path(datadir, "Figure3", "macrophage_module_scores.csv")
if (file.exists(score_path)) {
  scores <- read.csv(score_path, check.names = FALSE)
  p <- ggplot(scores, aes(x = Mac, y = score, fill = Mac, color = Mac)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.25) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_wrap(~score_name, scales = "free_y") +
    scale_fill_manual(values = mac_palette) +
    scale_color_manual(values = mac_palette) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = "", y = "Score")
  save_panel("Figure3D-G_macrophage_module_scores.pdf", p, width = 8, height = 6)
}

# Figure 3F/S3E: SCENIC/regulon activity plots
scenic_plot <- file.path(datadir, "Figure3", "macrophage_scenic_plot.rds")
if (file.exists(scenic_plot)) {
  p <- readRDS(scenic_plot)
  save_panel("Figure3F_macrophage_SCENIC.pdf", p, width = 6, height = 5)
}

t_path <- file.path(datadir, "Figure3", "T_cell.rds")
if (file.exists(t_path)) {
  T_cell <- readRDS(t_path)
  Idents(T_cell) <- "T"

  # Figure 3H: T-cell UMAP
  p <- DimPlot(T_cell, reduction = "umap", group.by = "T", cols = t_palette[seq_len(nlevels(T_cell))], label = TRUE, repel = TRUE, raster = TRUE) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  save_panel("Figure3H_Tcell_umap.pdf", p, width = 5, height = 4)

  # Figure 3K: type II IFN production score
  if ("type_ii_interferon_production1" %in% colnames(T_cell@meta.data)) {
    p <- VlnPlot(T_cell, features = "type_ii_interferon_production1", group.by = "T", pt.size = 0, cols = t_palette) +
      theme_classic() +
      labs(x = "", y = "Type II IFN production score")
    save_panel("Figure3K_Tcell_typeII_IFN_score.pdf", p, width = 6, height = 4)
  }
}

# Figure 3I: T subtype pathway enrichment barplot
t_enrich_path <- file.path(datadir, "Figure3", "Tcell_GO_terms.csv")
if (file.exists(t_enrich_path)) {
  go <- read.csv(t_enrich_path, check.names = FALSE)
  p <- ggplot(go, aes(x = reorder(Description, Count), y = Count, fill = subtype)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~subtype, scales = "free_y") +
    theme_classic() +
    labs(x = "", y = "Gene count")
  save_panel("Figure3I_Tcell_GO_barplot.pdf", p, width = 8, height = 6)
}

# Figure 3J: NicheNet ligands for CD8_Tef
nichenet_plot <- file.path(datadir, "Figure3", "CD8_Tef_nichenet_plot.rds")
if (file.exists(nichenet_plot)) {
  p <- readRDS(nichenet_plot)
  save_panel("Figure3J_CD8_Tef_NicheNet.pdf", p, width = 8, height = 4)
}
