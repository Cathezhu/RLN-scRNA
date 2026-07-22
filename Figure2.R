suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
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

time_levels <- c("Control", "D3", "D7", "D14", "D28")
time_palette <- c("Control" = "#1B9E77", "D3" = "#D95F02", "D7" = "#7570B3", "D14" = "#E7298A", "D28" = "#66A61E")

# Figure 2A-D: representative temporal gene trajectories
# Expected columns: gene, time, avg_exp, group
trajectory_path <- file.path(datadir, "Figure2", "temporal_gene_trajectories.csv")
if (file.exists(trajectory_path)) {
  traj <- read.csv(trajectory_path, check.names = FALSE)
  traj$time <- factor(traj$time, levels = time_levels)
  p <- ggplot(traj, aes(x = time, y = avg_exp, group = gene)) +
    geom_line(color = "#CD9963", linewidth = 1) +
    geom_point(color = "#CD9963", size = 2) +
    facet_wrap(~group + gene, scales = "free_y") +
    theme_classic() +
    labs(x = "", y = "Average expression")
  save_panel("Figure2A-D_temporal_gene_trajectories.pdf", p, width = 10, height = 8)
}

# Figure 2E: heatmap of time-point-specific genes
# Expected matrix rows: genes; columns: time or cell/time summaries.
heatmap_path <- file.path(datadir, "Figure2", "timepoint_gene_heatmap_matrix.csv")
if (file.exists(heatmap_path)) {
  mat <- read.csv(heatmap_path, row.names = 1, check.names = FALSE) %>% as.matrix()
  mat <- t(scale(t(mat)))
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  pdf(file.path(figdir, "Figure2E_timepoint_gene_heatmap.pdf"), width = 6, height = 8)
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000"))(100))
  dev.off()
}

# Figure 2G-H: GO/module scores across cell types/time points
# Expected columns: score_name, celltype, time, score
score_path <- file.path(datadir, "Figure2", "module_scores_long.csv")
if (file.exists(score_path)) {
  scores <- read.csv(score_path, check.names = FALSE)
  scores$time <- factor(scores$time, levels = time_levels)
  p <- ggplot(scores, aes(x = time, y = score, fill = time)) +
    geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_grid(score_name ~ celltype, scales = "free_y") +
    scale_fill_manual(values = time_palette) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = "", y = "Module score")
  save_panel("Figure2G-H_module_scores.pdf", p, width = 12, height = 8)
}

# Figure 2I: dynamic trend heatmap from Mfuzz/sub-Mfuzz outputs
trend_path <- file.path(datadir, "Figure2", "dynamic_trend_heatmap_matrix.csv")
if (file.exists(trend_path)) {
  mat <- read.csv(trend_path, row.names = 1, check.names = FALSE) %>% as.matrix()
  col_fun <- circlize::colorRamp2(seq(-2, 2, length = 100), colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100))
  pdf(file.path(figdir, "Figure2I_dynamic_trend_heatmap.pdf"), width = 6, height = 8)
  ComplexHeatmap::Heatmap(mat, name = "z-score", col = col_fun, cluster_rows = TRUE, cluster_columns = FALSE)
  dev.off()
}
