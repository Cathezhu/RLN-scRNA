suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

# ------------- Figure 4 --------------
# Plotting code for fibroblast/myofibroblast and predicted macrophage-fibroblast communication panels.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

fib_palette <- c("#3CB371", "#A97DBA", "#F89C1E", "#A2B7C7", "#7DDCFB", "#FB9A99", "#99003B", "#003366", "#B49C8F", "#9C9C9C", "#EC3A8D")

fib_path <- file.path(datadir, "Figure4", "Mesenchymal_cell.rds")
if (file.exists(fib_path)) {
  Fib <- readRDS(fib_path)
  Idents(Fib) <- "Fib"

  # Figure 4A: fibroblast/epineurial UMAP
  p <- DimPlot(Fib, reduction = "umap", group.by = "Fib", cols = fib_palette[seq_len(nlevels(Fib))], label = TRUE, repel = TRUE, raster = TRUE) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  save_panel("Figure4A_fibroblast_umap.pdf", p, width = 5, height = 4)

  # Figure 4B: fibroblast subtype composition
  meta <- Fib@meta.data %>%
    count(orig.ident, Fib, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()
  p <- ggplot(meta, aes(x = orig.ident, y = percent, fill = Fib)) +
    geom_col(position = "fill", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = fib_palette) +
    theme_classic() +
    labs(x = "", y = "Fraction")
  save_panel("Figure4B_fibroblast_fraction.pdf", p, width = 4, height = 4)

  # Figure 4C: fibroblast marker dotplot
  mark_feature <- list(
    Fib_1 = c("COL21A1", "SEMA3B", "SEMA3E", "VIT", "CPXM2"),
    Fib_2 = c("RSAD2", "IFIT3", "ISG15", "CXCL10"),
    Fib_3 = c("COL24A1", "TNC", "PTN", "TGFBI", "SFRP1"),
    Fib_4 = c("COL15A1", "ROBO2", "CCL19", "CXCL14", "IL6", "IGF1"),
    Fib_5 = c("ANXA8", "MFAP5", "PI16", "FBLN2"),
    Myofibroblast = c("ACTA2", "COL8A1", "POSTN", "COL27A1", "TGFB1")
  )
  mark_feature <- lapply(mark_feature, intersect, y = rownames(Fib))
  p <- DotPlot(Fib, features = mark_feature, group.by = "Fib") +
    coord_flip() +
    theme_bw() +
    labs(x = "", y = "")
  save_panel("Figure4C_fibroblast_marker_dotplot.pdf", p, width = 8.5, height = 3)
}

# Figure 4D: fibroblast GO dotplot
go_plot <- file.path(datadir, "Figure4", "fibroblast_GO_dotplot.rds")
if (file.exists(go_plot)) {
  p <- readRDS(go_plot)
  save_panel("Figure4D_fibroblast_GO.pdf", p, width = 8, height = 5)
}

# Figure 4E: PROGENy pathway activity
progeny_plot <- file.path(datadir, "Figure4", "fibroblast_PROGENy_plot.rds")
if (file.exists(progeny_plot)) {
  p <- readRDS(progeny_plot)
  save_panel("Figure4E_fibroblast_PROGENy.pdf", p, width = 7, height = 5)
}

# Figure 4F: Slingshot trajectory
slingshot_plot <- file.path(datadir, "Figure4", "fibroblast_slingshot_plot.rds")
if (file.exists(slingshot_plot)) {
  p <- readRDS(slingshot_plot)
  save_panel("Figure4F_fibroblast_slingshot.pdf", p, width = 6, height = 5)
}

# Figure 4G: collagen/ECM scores
ecm_score_path <- file.path(datadir, "Figure4", "ECM_scores_long.csv")
if (file.exists(ecm_score_path)) {
  scores <- read.csv(ecm_score_path, check.names = FALSE)
  p <- ggplot(scores, aes(x = group, y = score, fill = group, color = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.25) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_wrap(~score_name, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(x = "", y = "Score")
  save_panel("Figure4G_ECM_scores.pdf", p, width = 8, height = 5)
}

# Figure 4H-J: CellChat heatmap/rankNet/SPP1 bubble.
# These are plotting-only from precomputed CellChat plot objects or exported data.
for (panel in c("Figure4H_CellChat_heatmap", "Figure4I_rankNet", "Figure4J_SPP1_bubble")) {
  plot_path <- file.path(datadir, "Figure4", paste0(panel, ".rds"))
  if (file.exists(plot_path)) {
    p <- readRDS(plot_path)
    save_panel(paste0(panel, ".pdf"), p, width = 6, height = 5)
  }
}
