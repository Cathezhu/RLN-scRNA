suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
})

# ------------- Figure 2 --------------
# Plotting code for temporal gene dynamics and pathway-score panels.
# Required plotting inputs are precomputed tables under ./data/Figure2.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

####Figure2E####
# Figure 2E: time-specific gene heatmap
memb_df_path <- file.path(datadir, "Figure2", "memb_df.csv")
avg_all_path <- file.path(datadir, "Figure2", "avg_all.txt")
if (file.exists(memb_df_path) && file.exists(avg_all_path) && exists("expr_matrix")) {
  memb_df <- read.csv(memb_df_path) %>% filter(source %in% paste0("Group", c(1, 4, 3, 6))) %>%
    mutate(source = factor(source, levels = paste0("Group", c(1, 4, 3, 6)))) %>% arrange(source, desc(MEM.SHIP))
  gene.matrix <- read.table(avg_all_path)
  expr_matrix <- cbind(expr_matrix %>% as.data.frame(),
                       gene = gene.matrix$Gene.Symbol,
                       cell = gene.matrix$cell)
  expr_matrix_select <- expr_matrix[memb_df$NAME, ]
  rownames(memb_df) <- memb_df$X
  expr_matrix_select$group <- memb_df$source

  expr_matrix_select <- expr_matrix_select[!duplicated(expr_matrix_select$gene), ]
  rownames(expr_matrix_select) <- expr_matrix_select$gene

  mark_genes_path <- file.path(datadir, "Figure2", "Figure2E_mark_genes.txt")
  MarkGenes <- if (file.exists(mark_genes_path)) {
    trimws(readLines(mark_genes_path, warn = FALSE))
  } else {
    character(0)
  }
  MarkGenes <- MarkGenes[nzchar(MarkGenes)]

  markerGeneL <- rownames(expr_matrix_select) %in% MarkGenes
  markerGenes <- rownames(expr_matrix_select)[markerGeneL]
  rowHa <- NULL
  if (length(markerGenes) > 0) {
    rowHa <- rowAnnotation(
      link = anno_mark(
        at = which(markerGeneL),
        labels = markerGenes, labels_gp = gpar(fontsize = 6, fontface = "italic"), link_gp = gpar(lwd = 0.2), padding = 0.5
      ),
      width = unit(0.2, "cm") + max_text_width(markerGenes, gp = gpar(fontsize = 6))
    )
  }

  ht_opt$TITLE_PADDING <- unit(c(2.5, 2.5), "points")
  cols <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
  breaks <- seq(-2, 2, length = 100)
  col1 <- circlize::colorRamp2(breaks, cols)

  pdf(file.path(figdir, "Figure2E.pdf"), width = 4.5, height = 6)
  Heatmap(as.matrix(expr_matrix_select[1:5]),
          name = "z-score",
          right_annotation = rowHa,
          cluster_rows = F, cluster_columns = F,
          col = col1, border = F,
          show_row_names = FALSE, show_column_names = TRUE,
          use_raster = TRUE, raster_resize_mat = max
  )

  dev.off()
}

####Figure2I####
# Figure 2I: dynamic trend heatmap from Mfuzz/sub-Mfuzz outputs
trend_path <- file.path(datadir, "Figure2", "dynamic_trend_heatmap_matrix.csv")
if (file.exists(trend_path)) {
  mat <- read.csv(trend_path, row.names = 1, check.names = FALSE) %>% as.matrix()
  col_fun <- circlize::colorRamp2(seq(-2, 2, length = 100), colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100))
  pdf(file.path(figdir, "Figure2I.pdf"), width = 6, height = 8)
  ComplexHeatmap::Heatmap(mat, name = "z-score", col = col_fun, cluster_rows = TRUE, cluster_columns = FALSE)
  dev.off()
}

####FigS2A####
# Fig S2A: GO enrichment dotplot
go_path <- file.path(datadir, "Figure2", "go.xlsx")
if (file.exists(go_path)) {
  enrich_example <- readxl::read_excel(go_path)
  enrich_example <- enrich_example %>%
    filter(select == 1, Cluster %in% paste0("Group", c(1, 4, 3, 6))) %>%
    mutate(
      Description = paste0(
        toupper(substr(Description, 1, 1)),
        substr(Description, 2, nchar(Description))
      )
    )
  enrich_example$Cluster <- factor(enrich_example$Cluster, levels = paste0("Group", c(1, 4, 3, 6)))
  p <- scplotter::EnrichmentPlot(
    enrich_example,
    plot_type = "comparison",
    group_by = "Cluster",
    top_term = 100,
    aspect.ratio = 2.8
  )
  save_panel("FigS2A.pdf", p, width = 15, height = 15)
}

####FigS2B####
# Fig S2B: ribosomal gene percentage violin plot
s2b_time_levels <- c("intact", "3dpi", "7dpi", "14dpi", "28dpi")
s2b_time_palette <- c("#DC4E32", "#45B6C4", "#259C87", "#3A548C", "#EA987B")
names(s2b_time_palette) <- s2b_time_levels
s2b_path <- file.path(datadir, "Figure2", "FigS2B_percent_rb.csv")
s2b_df <- NULL

if (exists("new_celltype")) {
  if (!"percent.rb" %in% colnames(new_celltype@meta.data)) {
    new_celltype[["percent.rb"]] <- Seurat::PercentageFeatureSet(new_celltype, pattern = "^RP[SL]")
  }
  s2b_df <- new_celltype@meta.data %>%
    dplyr::select(time = orig.ident, percent.rb)
} else if (file.exists(s2b_path)) {
  s2b_df <- read.csv(s2b_path, check.names = FALSE)
}

if (!is.null(s2b_df)) {
  s2b_df$time <- factor(s2b_df$time, levels = s2b_time_levels)

  p <- ggplot(s2b_df, aes(x = time, y = percent.rb, fill = time)) +
    geom_violin(scale = "width", trim = FALSE, color = "black", linewidth = 0.35) +
    scale_fill_manual(values = s2b_time_palette) +
    scale_y_continuous(breaks = seq(0, 40, 10), expand = expansion(mult = c(0, 0.03))) +
    coord_cartesian(ylim = c(0, 40)) +
    labs(x = "", y = "", title = "percent.rb") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      legend.position = "none"
    )
  save_panel("FigS2B.pdf", p, width = 3.2, height = 2.6)
}
