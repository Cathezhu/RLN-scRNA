suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

# ------------- Figure 1 --------------
# Plotting code for global single-cell atlas panels.
# Required processed input:
#   ./data/Figure1/new_celltype.rds
# Optional inputs:
#   ./data/Figure1/markers_feature_heatmap.csv
#   ./data/Figure1/irGSEA_heatmap_plot.rds

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

global_palette <- c(
  "#E64B35FF", "#f4c7b4", "#F39B7FFF",
  "#fae6b8", "#dfad80", "#B09C85FF",
  "#91D1C2FF", "#00A087FF", "#8491B4FF",
  "#ae93b8", "#4DBBD5FF", "#3C5488FF",
  "#7E6148FF", "#98c98b"
)

load_obj <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  readRDS(path)
}

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

obj <- load_obj(file.path(datadir, "Figure1", "new_celltype.rds"))
Idents(obj) <- "new_celltype"
pal_use <- rep(global_palette, length.out = length(levels(obj)))

# Figure 1E: integrated UMAP of major cell types
p <- DimPlot(obj, reduction = "umap", group.by = "new_celltype", cols = pal_use, label = TRUE, repel = TRUE, raster = TRUE) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
save_panel("Figure1E_umap.pdf", p, width = 6, height = 5)

# Figure 1F: UMAP split by time point
p <- DimPlot(obj, reduction = "umap", group.by = "new_celltype", split.by = "orig.ident", cols = pal_use, raster = TRUE, ncol = 5) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
save_panel("Figure1F_umap_split.pdf", p, width = 12, height = 4)

# Figure 1G: cell-type composition across time points
meta <- obj@meta.data %>%
  count(orig.ident, new_celltype, name = "n") %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
meta$new_celltype <- factor(meta$new_celltype, levels = levels(obj))
p <- ggplot(meta, aes(x = orig.ident, y = new_celltype)) +
  geom_point(aes(size = percent, color = percent)) +
  scale_color_viridis_c(option = "magma", direction = 1) +
  scale_size(range = c(0.5, 7)) +
  theme_bw() +
  labs(x = "", y = "", color = "Percent", size = "Percent")
save_panel("Figure1G_composition_dot.pdf", p, width = 5, height = 5)

# Figure 1H: blue multi-gene UMAP feature plot
mark_feature <- c(
  "AIF1", "ITGAM", "CD68", "CD163", "C1QA", "MRC1",
  "CSF1R", "CD14", "S100A4", "FCGR3A", "FCN1",
  "LAMP3", "BCL11A", "FLT3", "FCER1A",
  "CD3E", "CD3D", "CD3G", "MS4A1", "CD79B", "CD19", "CD79A",
  "JCHAIN", "MZB1", "PDGFRA", "PDGFRB", "COL1A1",
  "BMP7", "SOX9", "ANGPTL7", "CDH5", "ENG", "PECAM1",
  "MYH11", "ACTA2", "RGS5", "MPZ", "MBP", "SOX10"
)
mark_feature <- intersect(mark_feature, rownames(obj))
umap_df <- Embeddings(obj, reduction = "umap") %>% as.data.frame()
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
plot_df <- data.frame(UMAP1 = umap_df$UMAP1, UMAP2 = umap_df$UMAP2, FetchData(obj, vars = mark_feature)) %>%
  pivot_longer(!c(UMAP1, UMAP2), names_to = "Markers", values_to = "Expr")
plot_df$Markers <- factor(plot_df$Markers, levels = mark_feature)
blue_cols <- RColorBrewer::brewer.pal(9, "Blues")
p <- ggplot(plot_df) +
  geom_point(aes(UMAP1, UMAP2, color = Expr), size = 0.5, stroke = 0, shape = 16) +
  facet_wrap(~Markers, scales = "free", ncol = 7) +
  theme_void() +
  theme(aspect.ratio = 1, strip.text = element_text(size = 8), legend.position = "bottom") +
  scale_color_gradientn(colours = blue_cols[c(1, 2, 5, 6, 7, 8, 9)], limits = c(0, 2.5), oob = scales::squish)
save_panel("Figure1H_marker_feature_blue.pdf", p, width = 10, height = 10, dpi = 500)

# Figure 1I: phenotype preference heatmap
# Store the finished ggplot/ComplexHeatmap object as ./data/Figure1/irGSEA_heatmap_plot.rds.
ir_plot_path <- file.path(datadir, "Figure1", "irGSEA_heatmap_plot.rds")
if (file.exists(ir_plot_path)) {
  p <- readRDS(ir_plot_path)
  save_panel("Figure1I_celltype_phenotype_heatmap.pdf", p, width = 10, height = 8)
}
