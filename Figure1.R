suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
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
new_celltype <- obj
Idents(new_celltype) <- "new_celltype"
celltype_palette <- rep(global_palette, length.out = nlevels(new_celltype))

####Figure1E####
# Figure 1E: integrated UMAP of major cell types
new_celltype <- subset(new_celltype, new_celltype %in% levels(new_celltype))
p <- scplotter::CellDimPlot(
  new_celltype,
  group_by = "new_celltype",
  reduction = "umap",
  label = TRUE,
  label_insitu = TRUE,
  theme = "theme_blank",
  palcolor = celltype_palette,
  label_size = 6,
  label_fg = "black",
  shuffle = TRUE,
  label_bg = "black",
  label_bg_r = 0,
  pt_size = 0.1,
  label_repel = TRUE,
  label_repulsion = 80
)
save_panel("Figure1E.pdf", p, width = 10, height = 6)
obj <- new_celltype

####Figure1F####
# Figure 1F: UMAP split by time point
p <- scplotter::CellDimPlot(
  new_celltype,
  group_by = "new_celltype",
  reduction = "umap",
  pt_size = 0.6,
  label = FALSE,
  split_by = "orig.ident",
  ncol = 2,
  theme = "theme_blank",
  combine = TRUE,
  legend.position = "none",
  palcolor = celltype_palette
)
save_panel("Figure1F.pdf", p, width = 12, height = 6)

####Figure1G####
# Figure 1G: cell-type composition across time points
p <- SCP::CellStatPlot(
  new_celltype,
  stat.by = "new_celltype",
  group.by = "orig.ident",
  label = TRUE,
  plot_type = "dot",
  palcolor = celltype_palette
)
save_panel("Figure1G.pdf", p, width = 5, height = 5)

####Figure1H####
# Figure 1H: rainbow marker UMAP feature plot
feature_palette <- RColorBrewer::brewer.pal(11, "Spectral") %>% rev()
mark_feature <- c(
  # Macrophage
  "CD68",
  # Monocyte
  "FCER1A",
  # T cell
  "CD3E",
  # B cell
  "CD79B",
  # Plasma cell
  "JCHAIN",
  # Mesenchymal cell
  "PDGFRA",
  # Endothelial cell
  "PECAM1",
  # Mural cell (SMC & pericyte)
  "MYH11",
  # Schwann cell
  "MPZ"
)
mark_feature <- intersect(mark_feature, rownames(obj))
umap_df <- Embeddings(obj, reduction = "umap") %>% as.data.frame()
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
plot_df <- data.frame(UMAP1 = umap_df$UMAP1, UMAP2 = umap_df$UMAP2, FetchData(obj, vars = mark_feature)) %>%
  pivot_longer(!c(UMAP1, UMAP2), names_to = "Markers", values_to = "Expr") %>%
  mutate(Expr = ifelse(Expr == 0, NA, Expr))
plot_df$Markers <- factor(plot_df$Markers, levels = mark_feature)
p <- ggplot(plot_df) +
  geom_point(aes(UMAP1, UMAP2, color = Expr), size = 0.8, stroke = 0, shape = 16) +
  facet_wrap(~Markers, scales = "free", ncol = 3) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    plot.margin = grid::unit(c(-2, -2, -2, -2), "char"),
    text = element_text(family = "sans"),
    strip.text = element_text(size = 8)
  ) +
  labs(color = "Exp") +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    legend.key.size = grid::unit(0.15, "inch"),
    legend.box.margin = margin(-10, 0, 0, 0)
  ) &
  scale_color_gradientn(
    colours = feature_palette,
    limits = c(0, 4),
    na.value = "grey90",
    oob = scales::squish
  )
save_panel("Figure1H.pdf", p, width = 5, height = 5, dpi = 500)

####Figure1I####
# Figure 1I: cell-type marker enrichment heatmap from precomputed irGSEA results
plot_irGSEA_heatmap <- function(object,
                                method = "RRA",
                                top = 50,
                                cluster_rows = FALSE,
                                cluster_color = NULL,
                                cluster_levels = NULL,
                                geneset_levels = NULL,
                                significance_color = c("#F4E1C1", "#ED4C4C"),
                                rowname_fontsize = 7,
                                heatmap_width = 10,
                                heatmap_height = 12) {
  if (!is.list(object)) {
    stop("object should be a list.")
  }
  if (!method %in% names(object)) {
    stop("method not found in object: ", method)
  }

  result <- object[[method]]
  if (!"pvalue" %in% colnames(result) && "p_val_adj" %in% colnames(result)) {
    result <- result %>% rename(pvalue = p_val_adj)
  }

  result <- result %>%
    filter(direction == "up") %>%
    mutate(
      Name = as.character(Name),
      cell = stringr::str_c(cluster, direction, sep = "_")
    )

  if (is.null(cluster_levels)) {
    cluster_levels <- unique(result$cluster)
  }
  if (is.null(cluster_color)) {
    cluster_color <- ggsci::pal_igv()(length(cluster_levels))
  }

  expected_cells <- paste0(cluster_levels, "_up")
  if (is.null(geneset_levels)) {
    geneset_levels <- unique(result$Name)
  }
  geneset_levels <- geneset_levels[geneset_levels %in% result$Name]
  geneset_levels <- head(geneset_levels, top)

  result <- result %>%
    filter(Name %in% geneset_levels) %>%
    mutate(Name = factor(Name, levels = geneset_levels))

  pvalue_table <- result %>%
    select(Name, pvalue, cell) %>%
    distinct() %>%
    tidyr::pivot_wider(names_from = cell, values_from = pvalue, values_fn = min) %>%
    arrange(Name) %>%
    as.data.frame()

  text_table <- result %>%
    mutate(label = case_when(
      pvalue < 1e-4 ~ "****",
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      is.na(pvalue) ~ "-",
      TRUE ~ "-"
    )) %>%
    select(Name, label, cell) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = cell,
      values_from = label,
      values_fill = "-",
      values_fn = list(label = function(x) x[which.max(nchar(x))])
    ) %>%
    arrange(Name) %>%
    as.data.frame()

  for (cell_name in setdiff(expected_cells, colnames(pvalue_table))) {
    pvalue_table[[cell_name]] <- NA_real_
  }
  for (cell_name in setdiff(expected_cells, colnames(text_table))) {
    text_table[[cell_name]] <- "-"
  }

  pvalue_table <- pvalue_table[, c("Name", expected_cells), drop = FALSE]
  text_table <- text_table[, c("Name", expected_cells), drop = FALSE]
  rownames(pvalue_table) <- as.character(pvalue_table$Name)
  rownames(text_table) <- as.character(text_table$Name)
  pvalue_table$Name <- NULL
  text_table$Name <- NULL

  heatmap_matrix <- -log10(as.matrix(pvalue_table))
  heatmap_matrix[is.na(heatmap_matrix)] <- 1
  heatmap_text <- as.matrix(text_table)

  significance_fun <- circlize::colorRamp2(
    c(1, 2, 3),
    colorRampPalette(significance_color)(3)
  )

  bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = ComplexHeatmap::anno_points(
      x = rep(1, length(expected_cells)),
      gp = grid::gpar(col = cluster_color, fill = cluster_color),
      pch = 16,
      size = grid::unit(3, "mm"),
      axis = FALSE,
      background_gp = grid::gpar(col = NA, fill = NA)
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 8),
    which = "column",
    annotation_height = grid::unit(2, "mm")
  )

  heatmap_body <- ComplexHeatmap::Heatmap(
    heatmap_matrix,
    na_col = "#f0f0f0",
    bottom_annotation = bottom_anno,
    heatmap_width = grid::unit(heatmap_width, "cm"),
    heatmap_height = grid::unit(heatmap_height, "cm"),
    name = method,
    col = significance_fun,
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    color_space = "RGB",
    show_column_names = FALSE,
    row_names_side = "right",
    row_names_max_width = ComplexHeatmap::max_text_width(
      rownames(heatmap_matrix),
      gp = grid::gpar(fontsize = rowname_fontsize)
    ),
    row_names_gp = grid::gpar(fontsize = rowname_fontsize),
    rect_gp = grid::gpar(col = "white", lwd = 2),
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(heatmap_text[i, j], x, y, gp = grid::gpar(fontsize = 10))
    }
  )

  pvalue_legend <- ComplexHeatmap::Legend(
    col_fun = significance_fun,
    title = "p-value",
    at = c(1, 2, 3, 4),
    labels = c("0.05", "0.01", "0.001", "0.0001"),
    legend_height = grid::unit(2, "cm")
  )
  significance_legend <- ComplexHeatmap::Legend(
    pch = c("-", "*", "**", "***", "****"),
    type = "points",
    labels = c("no significance", "< 0.05", "< 0.01", "< 0.001", "< 0.0001"),
    title = "P Value"
  )
  heatmap_legend <- ComplexHeatmap::packLegend(
    pvalue_legend,
    significance_legend,
    direction = "vertical",
    column_gap = grid::unit(1, "cm")
  )

  grid::grid.grabExpr(
    ComplexHeatmap::draw(heatmap_body, annotation_legend_list = heatmap_legend)
  ) %>%
    ggplotify::as.ggplot()
}

result_dge_path <- file.path(datadir, "Figure1", "result.dge.rds")
if (file.exists(result_dge_path)) {
  result_dge <- readRDS(result_dge_path)
  geneset_order <- unique(result_dge[[1]]$Name)
  if (length(geneset_order) >= 16) {
    geneset_order <- geneset_order[c(7, 4, 1, 5, 2, 3, 8:12, 6, 15, 16, 13, 14)]
  }
  p <- plot_irGSEA_heatmap(
    object = result_dge,
    method = "RRA",
    top = 50,
    cluster_rows = FALSE,
    cluster_color = celltype_palette,
    cluster_levels = levels(new_celltype),
    geneset_levels = geneset_order,
    significance_color = c("#F4E1C1", "#ED4C4C"),
    heatmap_width = 10,
    heatmap_height = 12
  ) +
    coord_flip()
  save_panel("Figure1I.pdf", p, width = 10, height = 8, dpi = 300)
}

####FigS1A####
# Fig S1A: unannotated global clusters
p <- DimPlot(
  new_celltype,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
save_panel("FigS1A.pdf", p, width = 5, height = 4)

####FigS1B-C####
# Fig S1B-C: QC metrics
qc_features <- intersect(c("nFeature_RNA", "percent.mt"), colnames(new_celltype@meta.data))
if (length(qc_features) > 0) {
  p <- VlnPlot(
    new_celltype,
    features = qc_features,
    group.by = "new_celltype",
    pt.size = 0,
    cols = celltype_palette,
    ncol = length(qc_features)
  ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_panel("FigS1B-C.pdf", p, width = 10, height = 4)
}

####FigS1D####
# Fig S1D: Zhang's blue multi-marker FeaturePlot
feature_palette <- RColorBrewer::brewer.pal(9, "Blues")
mark_feature <- c(
  "AIF1", "ITGAM",
  "CD68", "CD163", "C1QA", "MRC1",
  "CSF1R", "CD14", "S100A4", "FCGR3A", "FCN1",
  "LAMP3", "BCL11A", "FLT3", "FCER1A",
  "CD3E", "CD3D", "CD3G",
  "MS4A1", "CD79B", "CD19", "CD79A",
  "JCHAIN", "MZB1",
  "PDGFRA", "PDGFRB", "COL1A1",
  "BMP7", "SOX9", "ANGPTL7",
  "CDH5", "ENG", "PECAM1",
  "MYH11", "ACTA2", "RGS5",
  "MPZ", "MBP", "SOX10"
)
mark_feature <- intersect(mark_feature, rownames(new_celltype))
if (length(mark_feature) > 0) {
  umap_df <- Embeddings(new_celltype, reduction = "umap") %>% as.data.frame()
  plot_df <- data.frame(
    UMAP1 = umap_df$umap_1,
    UMAP2 = umap_df$umap_2,
    FetchData(object = new_celltype, vars = mark_feature)
  ) %>%
    tidyr::pivot_longer(
      cols = -c(UMAP1, UMAP2),
      names_to = "Markers",
      values_to = "Expr"
    )

  plot_df$Markers <- factor(plot_df$Markers, levels = mark_feature)

  p <- ggplot(plot_df) +
    geom_point(
      aes(x = UMAP1, y = UMAP2, color = Expr),
      size = 0.5,
      stroke = 0,
      shape = 16
    ) +
    theme_void() +
    theme(
      aspect.ratio = 1,
      plot.margin = grid::unit(c(-2, -2, -2, -2), "char"),
      text = element_text(family = "sans")
    ) +
    facet_wrap(~Markers, scales = "free", ncol = 7) +
    theme(strip.text = element_text(size = 8)) +
    labs(color = "Exp") +
    patchwork::plot_layout(guides = "collect") &
    theme(
      legend.title = element_text(size = 3.5),
      legend.text = element_text(size = 3),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.15, "inch"),
      legend.box.margin = margin(-10, 0, 0, 0)
    ) &
    scale_color_gradientn(
      colours = feature_palette[c(1, 2, 5, 6, 7, 8, 9)],
      limits = c(0, 2.5),
      oob = scales::squish
    )
  save_panel("FigS1D.pdf", p, width = 10, height = 10, dpi = 500)
}
