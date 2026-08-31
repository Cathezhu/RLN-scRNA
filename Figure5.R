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
    ####Figure5A####
    # Figure 5A: SPP1 expression across time
    p <- VlnPlot(Macrophage, features = "SPP1", group.by = "orig.ident", pt.size = 0) +
      theme_classic() +
      labs(x = "", y = "SPP1 expression")
    save_panel("Figure5A.pdf", p, width = 5, height = 4)

    ####Figure5B####
    # Figure 5B: SPP1 expression across macrophage subtype
    p <- VlnPlot(Macrophage, features = "SPP1", group.by = "Mac", pt.size = 0, cols = mac_palette) +
      theme_classic() +
      labs(x = "", y = "SPP1 expression")
    save_panel("Figure5B.pdf", p, width = 5, height = 4)
  }
}

####Figure5C####
# Figure 5C: D3 CellChat heatmap from precomputed plot
plot_path <- file.path(datadir, "Figure5", "Figure5C_D3_CellChat_heatmap.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("Figure5C.pdf", p, width = 6, height = 5)
}

####Figure5D####
# Figure 5D: SuperCell/metacell correlation
if (!exists("Macrophage")) {
  mac_supercell_path <- Sys.getenv(
    "RLN_MACROPHAGE_RDS",
    unset = file.path(datadir, "Figure5", "Macrophage.rds")
  )
  if (!file.exists(mac_supercell_path)) {
    mac_supercell_path <- "E:/10X_RLN/new_qc/sub/Macrophage/Macrophage1.rds"
  }
  if (file.exists(mac_supercell_path)) {
    Macrophage <- readRDS(mac_supercell_path)
  }
}
if (!exists("Macrophage")) {
  stop("Provide Macrophage or set RLN_MACROPHAGE_RDS for Figure5D.")
}
if (!requireNamespace("SuperCell", quietly = TRUE)) {
  stop("Figure5D requires the SuperCell package.")
}
if (!"package:SuperCell" %in% search()) {
  suppressPackageStartupMessages(library(SuperCell))
}

supercell_helper <- Sys.getenv("RLN_SUPERCELL_HELPER", unset = "supercell_2_Seurat.R")
if (!file.exists(supercell_helper)) {
  stop("Missing SuperCell Seurat helper: ", supercell_helper)
}
source(supercell_helper)

obj <- Macrophage
exp_mat <- GetAssayData(obj, layer = "data", assay = "RNA")
meta <- obj$orig.ident
obj <- FindVariableFeatures(obj, nfeatures = 7000)
hvg <- VariableFeatures(obj)
SC <- SCimplify(
  exp_mat,
  k.knn = 5,
  gamma = 10,
  n.var.genes = 1000
)

SC.GE <- supercell_GE(exp_mat, SC$membership)
SC$cell_line <- supercell_assign(
  clusters = meta,
  supercell_membership = SC$membership,
  method = "jaccard"
)
supercell_plot(
  SC$graph.supercells,
  group = SC$cell_line,
  seed = 1,
  main = "Metacells colored by cell line assignment"
)
purity <- supercell_purity(
  clusters = meta,
  supercell_membership = SC$membership,
  method = "entropy"
)
hist(purity, main = "Purity of metacells \nin terms of cell line composition")
SC$purity <- purity

MC.seurat <- supercell_2_Seurat(
  SC.GE = SC.GE,
  SC = SC,
  fields = c("cell_line", "purity"),
  var.genes = SC$genes.use,
  N.comp = 10
)
p <- scplotter::FeatureStatPlot(
  MC.seurat,
  features = c("SPP1", "FN1", "FABP4"),
  plot_type = "cor"
)
save_panel("Figure5D.pdf", p, width = 5, height = 4)

####Figure5E####
# Figure 5E: macrophage monocle pseudotime
if (!requireNamespace("monocle", quietly = TRUE)) {
  stop("Figure5E requires the monocle package.")
}
if (!exists("cds")) {
  cds_path <- Sys.getenv(
    "RLN_MACROPHAGE_MONOCLE_CDS",
    unset = file.path(datadir, "Figure5", "cds.rds")
  )
  if (!file.exists(cds_path)) {
    cds_path <- "E:/10X_RLN/new_qc/sub/Macrophage/monocle/cds.rds"
  }
  if (!file.exists(cds_path)) {
    stop("Provide cds or set RLN_MACROPHAGE_MONOCLE_CDS for Figure5E.")
  }
  cds <- readRDS(cds_path)
}

p <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.8) +
  scale_color_gradientn(
    colors = rev(c("#FCFFA4FF", "#F3E45CFF", "#F9A03F", "#85216BFF", "#230C4BFF")),
    na.value = "grey80"
  )
save_panel("Figure5E.pdf", p, width = 6, height = 4)

####Figure5F####
# Figure 5F: branch heatmap
if (!requireNamespace("ClusterGVis", quietly = TRUE)) {
  stop("Figure5F requires the ClusterGVis package.")
}

beam_path <- Sys.getenv(
  "RLN_MACROPHAGE_BEAM_RDS",
  unset = file.path(datadir, "Figure5", "BEAM_res_branch_1.rds")
)
if (!file.exists(beam_path)) {
  beam_path <- "E:/10X_RLN/new_qc/sub/Macrophage/monocle/BEAM_res_branch_1.rds"
}
df_path <- Sys.getenv(
  "RLN_MACROPHAGE_BEAM_DF",
  unset = file.path(datadir, "Figure5", "moncoleDF.rds")
)
if (!file.exists(df_path)) {
  df_path <- "E:/10X_RLN/new_qc/sub/Macrophage/monocle/moncoleDF.rds"
}
if (!file.exists(beam_path) || !file.exists(df_path)) {
  stop("Provide BEAM_res_branch_1.rds and moncoleDF.rds for Figure5F.")
}

BEAM_res <- readRDS(beam_path) %>%
  dplyr::arrange(qval) %>%
  dplyr::select(gene_short_name, pval, qval) %>%
  dplyr::filter(!grepl("^RP[SL]", gene_short_name))
diff_gene <- BEAM_res %>%
  dplyr::filter(qval < 1e-4) %>%
  dplyr::arrange(qval)
df <- readRDS(df_path)

grDevices::pdf(file.path(figdir, "Figure5F.pdf"), width = 7, height = 9)
ClusterGVis::visCluster(
  object = df,
  plot.type = "heatmap",
  pseudotime_col = c("#7A97C3", "grey50", "#D66B65"),
  markGenes = unique(head(diff_gene$gene_short_name, 100)),
  ht.col.list = list(col_range = c(-3, 0, 3)),
  genes.gp = c("italic", 4, NA),
  ctAnno.col = mac_palette,
  show_row_names = FALSE,
  line.col = "white",
  cluster.order = c(4, 5, 2, 1, 3)
)
grDevices::dev.off()

####Figure5G####
# Figure 5G: pseudotime gene dynamics
plot_path <- file.path(datadir, "Figure5", "Figure5G_pseudotime_gene_dynamics.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("Figure5G.pdf", p, width = 6, height = 5)
}

####Figure5H####
# Figure 5H: Mac_4 vs Mac_5 GSEA comparison
gsea_path <- Sys.getenv(
  "RLN_MAC4_VS_MAC5_GSEA",
  unset = file.path(datadir, "Figure5", "Mac4vs5.csv")
)
if (!file.exists(gsea_path)) {
  stop("Provide Mac4vs5 GSEA results or set RLN_MAC4_VS_MAC5_GSEA for Figure5H.")
}
GSEA_term_Sig <- read.csv(gsea_path, check.names = FALSE) %>%
  dplyr::filter(select == 1)
GSEA_term_Sig$Description <- GSEA_term_Sig$Description %>%
  stringr::str_remove("_GO:\\d+") %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_to_title()
GSEA_term_Sig$Logqvalue <- ifelse(
  GSEA_term_Sig$NES > 0,
  -log(GSEA_term_Sig$qvalue),
  log(GSEA_term_Sig$qvalue)
)
GSEA_term_Sig$Logqvalue <- Seurat::MinMax(GSEA_term_Sig$Logqvalue, min = -10, max = 10)
GSEA_term_Sig$Cluster <- ifelse(GSEA_term_Sig$NES > 0, "Mac_4", "Mac_5")
col1 <- mac_palette

p <- ggpubr::ggbarplot(
  GSEA_term_Sig,
  x = "Description",
  y = "Logqvalue",
  fill = "Cluster",
  color = "white",
  palette = col1,
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "-Log10(qvalue)",
  xlab = FALSE,
  width = 0.98,
  title = "SPP+ FN1+ Mac vs SPP+ FABP4+ Mac"
) +
  coord_flip() +
  xlab("") +
  theme(
    axis.text.y = element_text(size = 17),
    plot.title = element_text(size = 30),
    legend.position = "none"
  )
save_panel("Figure5H.pdf", p, width = 12, height = 6)

####FigS5F####
# Fig S5F: macrophage Slingshot lineages
if (!exists("Macrophage")) {
  mac_slingshot_path <- Sys.getenv(
    "RLN_MACROPHAGE_RDS",
    unset = file.path(datadir, "Figure5", "Macrophage.rds")
  )
  if (!file.exists(mac_slingshot_path)) {
    mac_slingshot_path <- "E:/10X_RLN/new_qc/sub/Macrophage/Macrophage1.rds"
  }
  if (file.exists(mac_slingshot_path)) {
    Macrophage <- readRDS(mac_slingshot_path)
  }
}
if (!exists("Macrophage")) {
  stop("Provide Macrophage or set RLN_MACROPHAGE_RDS for FigS5F.")
}
if (!requireNamespace("scop", quietly = TRUE)) {
  stop("FigS5F requires the scop package.")
}

pancreas_sub <- scop::RunSlingshot(
  srt = Macrophage,
  group.by = "RNA_snn_res.1",
  start = "10",
  reduction = "umap",
  extend = "y",
  shrink = FALSE
)
p <- scop::CellDimPlot(
  pancreas_sub,
  group.by = "Mac",
  reduction = "umap",
  lineages = paste0("Lineage", c(1, 3, 4)),
  lineages_span = 0.5,
  palcolor = c(
    "#4eb897", "#f08033", "#9993cb", "#f062ae",
    "#ddc9e5", "#fdb94c", "#a1c8e4", "#aee2d9"
  ),
  lineages_palette = "Set1",
  theme_use = "theme_blank",
  label = TRUE,
  label_insitu = TRUE
)
save_panel("FigS5F.pdf", p, width = 8, height = 6)

####FigS5G####
# Fig S5G: macrophage Slingshot dynamic genes
if (!exists("pancreas_sub")) {
  stop("FigS5G requires pancreas_sub created by FigS5F.")
}

p <- scop::DynamicPlot(
  pancreas_sub,
  lineages = c("Lineage4", "Lineage3"),
  group.by = "Mac",
  features = c("SPP1", "FABP4", "FN1", "CCL24"),
  compare_lineages = TRUE,
  nrow = 2,
  compare_features = FALSE,
  point_palcolor = c(
    "#4eb897", "#f08033", "#9993cb", "#f062ae",
    "#ddc9e5", "#fdb94c", "#a1c8e4", "#aee2d9"
  ),
  pt.size = 0.2
)
save_panel("FigS5G.pdf", p, width = 7, height = 5)

####FigS5H####
# Fig S5H: macrophage Slingshot panel
plot_path <- file.path(datadir, "Figure5", "FigS5H.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("FigS5H.pdf", p, width = 7, height = 5)
}
