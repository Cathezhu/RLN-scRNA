suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

# ------------- Figure 4 --------------
# Plotting code for fibroblast/myofibroblast and predicted macrophage-fibroblast communication panels.

datadir <- Sys.getenv("RLN_FIGURE_DATA", unset = "./data")
figdir <- Sys.getenv("RLN_FIGURE_OUTPUT", unset = "./figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(filename, plot, width, height, dpi = 300) {
  ggsave(file.path(figdir, filename), plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

yy_Dotplot <- function(seuratObj,
                       genes,
                       group.by,
                       coord_flip = FALSE,
                       scale = TRUE,
                       dot.scale = 4,
                       gene_expr_cutoff = -Inf,
                       gene_pct_cutoff = 0,
                       cell_expr_cutoff = -Inf,
                       cell_pct_cutoff = 0,
                       panel.spacing_distance = 0.5,
                       return_data = FALSE) {
  col.min <- -2.5
  col.max <- 2.5
  dot.min <- 0

  if (is.list(genes)) {
    if (is.null(names(genes))) {
      names(genes) <- seq_along(genes)
    }
    genes <- stack(genes)
    features <- genes$values
    features_group <- genes$ind
  } else {
    features <- genes
    features_group <- NULL
  }

  subset <- features %in% rownames(seuratObj)
  if (sum(!subset) > 0) {
    cat(paste0(paste(features[!subset], collapse = ", "), " is missing in gene list"))
  }
  if (length(features) == 0) {
    stop("No intersecting genes, please check gene name format.\n")
  }

  data.features <- FetchData(seuratObj, cells = colnames(seuratObj), vars = features, slot = "data")
  data.features$id <- seuratObj@meta.data[[group.by]]

  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, seq_len(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, function(x) sum(x > 0) / length(x))
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  names(data.plot) <- unique(data.features$id)
  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$features.plot_show <- features
    data.use$id <- x
    data.use
  })
  data.plot <- do.call(rbind, data.plot)

  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    if (scale) {
      data.use <- scale(data.use)
      Seurat::MinMax(data.use, min = col.min, max = col.max)
    } else {
      log1p(data.use)
    }
  })
  data.plot$avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  data.plot$features.plot <- factor(data.plot$features.plot, levels = unique(data.plot$features.plot))
  data.plot$features.plot_show <- factor(data.plot$features.plot_show, levels = unique(features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100

  if (is.factor(seuratObj@meta.data[[group.by]])) {
    data.plot$id <- factor(data.plot$id, levels = levels(seuratObj@meta.data[[group.by]]))
  }
  if (!is.null(features_group)) {
    data.plot <- data.plot %>%
      left_join(data.frame(features = features, features_group = features_group), by = c("features.plot_show" = "features")) %>%
      mutate(features_group = factor(features_group, levels = unique(features_group)))
  }

  filter_cell <- data.plot %>%
    group_by(id) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp), .groups = "drop") %>%
    filter(max_exp >= cell_expr_cutoff & max_pct >= cell_pct_cutoff)
  filter_gene <- data.plot %>%
    group_by(features.plot) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp), .groups = "drop") %>%
    filter(max_exp >= gene_expr_cutoff & max_pct >= gene_pct_cutoff)
  data.plot <- data.plot %>%
    filter(id %in% filter_cell$id, features.plot %in% filter_gene$features.plot)

  plot <- ggplot(data.plot, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
    guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.2)
    ) +
    labs(x = "", y = "")

  if (coord_flip) {
    plot <- plot + coord_flip()
  }
  if (!is.null(features_group)) {
    plot <- plot +
      facet_grid(. ~ features_group, scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(panel.spacing_distance, "lines"), strip.background = element_blank())
  }
  if (return_data) {
    return(list(plot = plot, data = data.plot))
  }
  plot
}

fib_palette <- c("#3CB371", "#A97DBA", "#F89C1E", "#A2B7C7", "#7DDCFB", "#FB9A99", "#99003B", "#003366", "#B49C8F", "#9C9C9C", "#EC3A8D")

fib_path <- file.path(datadir, "Figure4", "Mesenchymal_cell.rds")
if (file.exists(fib_path)) {
  Fib <- readRDS(fib_path)
  Idents(Fib) <- "Fib"

  ####Figure4A####
  # Figure 4A: fibroblast/epineurial UMAP
  p <- scplotter::CellDimPlot(
    Fib,
    group_by = "Fib",
    reduction = "umap",
    label = TRUE,
    label_insitu = TRUE,
    theme = "theme_blank",
    palcolor = fib_palette[seq_len(nlevels(Fib))],
    label_size = 6,
    label_fg = "black",
    shuffle = TRUE,
    label_bg = "black",
    label_bg_r = 0,
    pt_size = 0.1,
    label_repel = TRUE,
    label_repulsion = 80
  )
  save_panel("Figure4A.pdf", p, width = 5, height = 4)

  ####Figure4B####
  # Figure 4B: fibroblast marker dotplot
  mark_feature <- list(
    Fib_1 = c("col21a1", "sema3b", "sema3e", "vit", "cpxm2"),
    Fib_2 = c("rsad2", "ifit3", "isg15", "cxcl10"),
    Fib_3 = c("col24a1", "tnc", "ptn", "tgfbi", "sfrp1"),
    Fib_4 = c("col15a1", "robo2", "ccl19", "cxcl14", "il6", "igf1"),
    Fib_5 = c("anxa8", "mfap5", "pi16", "fbln2"),
    Myofibroblast = c("acta2", "col8a1", "postn", "col27a1", "tgfb1")
  ) %>%
    lapply(toupper)
  names(mark_feature) <- NULL

  p <- yy_Dotplot(
    seuratObj = Fib,
    genes = mark_feature,
    group.by = "Fib"
  )
  p <- p +
    theme(
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank(),
      text = element_text(size = 0),
      panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.2),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.key.size = unit(0.15, "inch"),
      plot.margin = unit(c(0, 0, 0, 0), "char")
    ) +
    labs(x = "", y = "") +
    scale_radius(limits = c(0, 100), range = c(0, 2.5))
  save_panel("Figure4B.pdf", p, width = 8.5, height = 2.7, dpi = 800)

  ####Figure4C####
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
  save_panel("Figure4C.pdf", p, width = 8.5, height = 3)
}

####Figure4D####
# Figure 4D: fibroblast GO enrichment barplot
if (exists("Fib")) {
  fib_go_path <- Sys.getenv(
    "RLN_FIB_GO",
    unset = file.path(datadir, "Figure4", "Hs_DEG_GO.csv")
  )
  if (!file.exists(fib_go_path)) {
    fib_go_path <- "E:/10X_RLN/new_qc/sub/non.immune/Mesenchymal_cell/markers/Hs_DEG_GO.csv"
  }
  if (!file.exists(fib_go_path)) {
    stop("Missing GO enrichment table for Figure4D: ", fib_go_path)
  }
  filtered_go <- read.csv(fib_go_path) %>% filter(select == 1)
  filtered_go$cluster <- factor(filtered_go$cluster, levels = levels(Fib))
  filtered_go$Description <- factor(filtered_go$Description, levels = unique(filtered_go$Description))
  fib_palette_use <- fib_palette[seq_len(nlevels(Fib))]
  names(fib_palette_use) <- levels(Fib)

  p <- ggplot(filtered_go, aes(x = Description, y = -log(qvalue), fill = cluster)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = fib_palette_use) +
    facet_grid(~cluster, scales = "free", space = "free_y") +
    theme_bw() +
    xlab("") +
    ylab("-log10(qvalue)") +
    theme(
      strip.text.x = element_text(colour = "black", angle = 360, size = 20, hjust = 0),
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 12),
      axis.ticks = element_line(linewidth = 1, colour = "black"),
      axis.line = element_line(linewidth = 1, colour = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_flip()
  save_panel("Figure4D.pdf", p, width = 15, height = 8)
}

####Figure4E####
# Figure 4E: PROGENy pathway activity
if (exists("Fib")) {
  options(RCurlOptions = list(ftp.use.epsv = FALSE))
  progeny_net_path <- Sys.getenv(
    "RLN_PROGENY_NET",
    unset = file.path(datadir, "Figure4", "progeny_net.rds")
  )
  if (!file.exists(progeny_net_path)) {
    progeny_net_path <- "E:/10X_RLN/new_qc/sub/T_cell/progeny_net.rds"
  }
  if (!file.exists(progeny_net_path)) {
    stop("Missing PROGENy network for Figure4E: ", progeny_net_path)
  }

  net <- readRDS(progeny_net_path)
  mat <- as.matrix(Fib[["Homo_RNA"]]$data)
  acts <- decoupleR::run_mlm(
    mat = mat,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = 5
  )
  Fib[["pathwaysmlm"]] <- acts %>%
    pivot_wider(
      id_cols = "source",
      names_from = "condition",
      values_from = "score"
    ) %>%
    tibble::column_to_rownames("source") %>%
    Seurat::CreateAssayObject(.)

  DefaultAssay(Fib) <- "pathwaysmlm"
  Fib <- ScaleData(Fib)
  Fib[["pathwaysmlm"]]$data <- Fib[["pathwaysmlm"]]$scale.data

  progeny_scores_df <- as.data.frame(t(GetAssayData(
    Fib,
    slot = "scale.data",
    assay = "pathwaysmlm"
  ))) %>%
    tibble::rownames_to_column("Cell") %>%
    pivot_longer(cols = -Cell, names_to = "Pathway", values_to = "Activity")
  CellsClusters <- data.frame(
    Cell = colnames(Fib),
    CellType = as.character(Fib@meta.data$Fib),
    stringsAsFactors = FALSE
  )
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

  summarized_progeny_scores <- progeny_scores_df %>%
    dplyr::group_by(Pathway, CellType) %>%
    dplyr::summarise(avg = mean(Activity), std = sd(Activity), .groups = "drop")
  sub_progeny_scores <- summarized_progeny_scores %>%
    filter(Pathway %in% c("TNFa", "NFkB", "JNK-STAT", "TGFb", "WNT"))

  cluster.lineages <- levels(Fib)
  sub_progeny_scores <- sub_progeny_scores %>%
    filter(CellType %in% cluster.lineages)
  sub_progeny_scores$CellType <- factor(sub_progeny_scores$CellType, levels = rev(cluster.lineages))
  sub_progeny_scores$Pathway <- factor(
    sub_progeny_scores$Pathway,
    levels = c("TNFa", "NFkB", "JNK-STAT", "TGFb", "WNT")
  )
  sub_progeny_scores <- sub_progeny_scores %>% arrange(Pathway)
  sub_progeny_scores$group <- sub_progeny_scores$CellType

  groupMean <- sub_progeny_scores %>%
    dplyr::group_by(Pathway) %>%
    dplyr::mutate(MeanV = mean(avg)) %>%
    dplyr::distinct(Pathway, .keep_all = TRUE) %>%
    dplyr::select(Pathway, MeanV)
  sub_progeny_scores$CellType <- factor(sub_progeny_scores$CellType)
  line_df <- sub_progeny_scores %>%
    distinct(group, Pathway, CellType) %>%
    mutate(y1 = -Inf, y2 = Inf)

  fib_palette_use <- fib_palette[seq_len(nlevels(Fib))]
  names(fib_palette_use) <- levels(Fib)
  p <- ggplot(sub_progeny_scores, aes(x = CellType, y = avg, colour = group)) +
    geom_segment(
      data = line_df,
      aes(x = CellType, xend = CellType, y = y1, yend = y2, colour = group),
      linetype = "solid",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(size = 2.4) +
    geom_hline(data = groupMean, aes(yintercept = MeanV), linetype = "dashed") +
    coord_flip() +
    scale_colour_manual(values = fib_palette_use) +
    facet_grid(group ~ Pathway, scales = "free") +
    theme_gray(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7),
      strip.background = element_blank()
    )
  save_panel("Figure4E.pdf", p, width = 6, height = 4)
}

####Figure4F####
# Figure 4F: Slingshot trajectory
if (exists("Fib")) {
  cell_pal <- function(cell_vars, pal_fun, ...) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100, ...)
      return(pal[cut(cell_vars, breaks = 100)])
    }

    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    pal[cell_vars]
  }

  sce <- Seurat::as.SingleCellExperiment(Fib, assay = "RNA")
  sce_slingshot1 <- slingshot::slingshot(
    sce,
    reducedDim = "UMAP",
    clusterLabels = "Fib",
    start.clus = "Fib_1",
    approx_points = 150,
    allow.breaks = FALSE,
    extend = "n",
    shrink = FALSE
  )
  sds <- slingshot::SlingshotDataSet(sce_slingshot1)
  cell_colors <- cell_pal(sce_slingshot1$Fib, scales::hue_pal())

  umap_df <- Embeddings(Fib, reduction = "umap") %>% as.data.frame()
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
  celltype_label <- cbind(umap_df, celltype = Fib$Fib) %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      UMAP1 = median(UMAP1),
      UMAP2 = median(UMAP2),
      .groups = "drop"
    )

  grDevices::pdf(file.path(figdir, "Figure4F.pdf"), width = 8, height = 6)
  plot(
    SingleCellExperiment::reducedDims(sce_slingshot1)[["UMAP"]],
    col = cell_colors,
    pch = 16,
    asp = 1,
    cex = 0.4
  )
  for (curve in sds@curves) {
    lines(curve$s, lwd = 2, col = "black", xpd = FALSE, ljoin = "bevel")
  }
  lines(sds, lwd = 2, col = "black", type = "lineage")
  text(
    x = celltype_label$UMAP1 - 1,
    y = celltype_label$UMAP2,
    labels = celltype_label$celltype
  )
  grDevices::dev.off()

  pancreas_sub <- scop::RunSlingshot(
    srt = Fib,
    group.by = "Fib",
    start = "Fib_1",
    stretch = 2,
    reduction = "umap",
    extend = "pc1",
    shrink = FALSE
  )
  pancreas_sub$Lineage1_3 <- dplyr::case_when(
    !is.na(pancreas_sub$Lineage3) ~ pancreas_sub$Lineage3,
    !is.na(pancreas_sub$Lineage1) ~ pancreas_sub$Lineage1
  )
  p <- scop::FeatureDimPlot(
    pancreas_sub,
    features = "Lineage3",
    reduction = "UMAP",
    lineages = "Lineage3",
    theme_use = "theme_blank",
    lineages_span = 0.5,
    lineages_palette = "Set1",
    combine = TRUE
  )
}

####Figure4G####
# Figure 4G: collagen/ECM glycoprotein and ECM regulator scores
if (!exists("Macrophage")) {
  mac_path <- Sys.getenv(
    "RLN_MACROPHAGE_RDS",
    unset = file.path(datadir, "Figure4", "Macrophage.rds")
  )
  if (!file.exists(mac_path)) {
    mac_path <- "E:/10X_RLN/new_qc/sub/Macrophage/Macrophage1.rds"
  }
  if (file.exists(mac_path)) {
    Macrophage <- readRDS(mac_path)
  }
}

if (exists("Macrophage") && exists("Fib")) {
  collagen_genes <- unique(toupper(c(
    "col31A2", "COL5A1", "COL5A2", "col31A1", "COL7A1", "col34A1", "COL3A1",
    "COL2A1", "COL6A1", "col31A2", "COL6A3", "COL6A2", "col34A1", "COL7A1", "COL5A1", "COL22A1",
    "col31A1", "COL4A3", "col32A1", "COL2A1", "COL3A1", "COL4A1", "COL27A1",
    "COL28A1", "COL6A5", "COL4A6", "col38A1", "COL25A1", "COL23A1", "COL4A5",
    "COL21A1", "COL9A1", "col36A1", "COL26A1", "COL6A6", "col33A1", "COL20A1",
    "col37A1", "COL24A1",
    "FN1", "FNDC1", "LAMA4", "MATN3", "NPNT", "FBLN2", "LAMA5", "MATN2", "NTNG1",
    "EYS", "PAPLN", "SNED1", "SLIT2", "NELL2", "THBS3", "RSPO3", "LAMB4", "IGFBP4",
    "ZP2", "CDCP2", "TINAG", "NELL1", "FBLN1", "NTNG2", "NTN5", "LGI3", "VWA1",
    "VWA7", "KCP", "POSTN", "VIT", "LAMC2", "NID1", "FBN2", "EMILIN1", "AGRN"
  )))
  ecm_regulator_genes <- toupper(c(
    "TGM2", "ADAM15", "ADAM23", "ITIH5", "PLOD1", "ADAMTS20", "ADAM22",
    "ADAMTSL4", "ADAMTS9", "ITIH3", "ADAMTS4", "CD109", "PAMR1", "ADAMTS10",
    "NGLY1", "PZP", "ADAMTS6", "EGLN1", "ADAMTSL3", "HPSE2", "P4HA2", "ITIH1",
    "PLAU", "P4HA3", "SERPINB5", "ADAMTSL1", "ADAM7", "LOXL3", "MMP16", "HABP2",
    "ADAM18", "P4HTM", "EGLN2", "OGFOD1", "HTRA3", "PCSK6", "MMP19", "PCSK5",
    "SERPINB7", "SPAM1", "PAPPA2", "MMP2", "SERPINB3", "C17orf58", "ADAMTS2",
    "A2ML1", "ADAMTS17"
  ))
  pathway_all <- list(
    "Collagen and ECM glycoprotein" = collagen_genes,
    "ECM regulator" = ecm_regulator_genes
  )
  pathway_order <- names(pathway_all)
  celltype_order <- c(paste0("Mac_", 1:6), paste0("Fib_", 1:5), "Myofibroblast")

  score_pathways <- function(object, group_col) {
    df <- object@meta.data[, c("orig.ident", group_col), drop = FALSE]
    names(df)[names(df) == group_col] <- "celltype"

    for (item in pathway_order) {
      features <- intersect(pathway_all[[item]], rownames(object))
      missing_features <- setdiff(pathway_all[[item]], features)
      if (length(missing_features) > 0) {
        message(item, " missing in ", deparse(substitute(object)), ": ", paste(missing_features, collapse = ", "))
      }
      if (length(features) == 0) {
        stop("No genes from ", item, " are available in ", deparse(substitute(object)), ".")
      }

      scored_object <- UCell::AddModuleScore_UCell(
        object,
        features = list(features),
        name = "addscore",
        assay = "RNA"
      )
      df[[item]] <- scored_object@meta.data[[ncol(scored_object@meta.data)]]
    }
    df
  }

  df_long <- dplyr::bind_rows(
    score_pathways(Macrophage, "Mac"),
    score_pathways(Fib, "Fib")
  ) %>%
    dplyr::filter(celltype %in% celltype_order) %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pathway_order),
      names_to = "sig",
      values_to = "score"
    )
  df_long$celltype <- factor(df_long$celltype, levels = celltype_order)
  df_long$sig <- factor(df_long$sig, levels = pathway_order)

  figure4g_palette <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#cab2d6", "#fda111",
    fib_palette[seq_len(6)]
  )
  names(figure4g_palette) <- celltype_order

  sig_mean <- df_long %>%
    dplyr::group_by(sig) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
  stat_df <- df_long %>%
    dplyr::group_by(celltype, sig) %>%
    dplyr::summarise(
      Q1 = quantile(score, 0.25, na.rm = TRUE),
      median = median(score, na.rm = TRUE),
      Q3 = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(df_long, aes(x = celltype, y = score, fill = celltype, color = celltype)) +
    geom_violin(
      scale = "width",
      adjust = 1,
      trim = FALSE
    ) +
    geom_hline(
      data = sig_mean,
      aes(yintercept = mean_score),
      linetype = "dashed",
      color = "grey25",
      linewidth = 0.4
    ) +
    geom_errorbar(
      data = stat_df,
      aes(x = celltype, ymin = Q1, ymax = Q3),
      width = 0.1,
      color = "black",
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    geom_segment(
      data = stat_df,
      inherit.aes = FALSE,
      aes(
        x = as.numeric(celltype) - 0.05,
        xend = as.numeric(celltype) + 0.05,
        y = median,
        yend = median
      ),
      color = "black",
      linewidth = 0.4
    ) +
    facet_wrap(~sig, ncol = 1, scales = "free_y", strip.position = "top", dir = "h") +
    theme_classic() +
    ylab("Gene score") +
    xlab("") +
    scale_fill_manual(values = figure4g_palette) +
    scale_color_manual(values = figure4g_palette) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      legend.position = "none",
      strip.text = element_text(size = 6),
      strip.background = element_blank()
    )
  save_panel("Figure4G.pdf", p, width = 10, height = 5)
}

####Figure4H####
# Figure 4H: D28 CellChat heatmap
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Figure4H requires the CellChat package.")
}
if (!"package:CellChat" %in% search()) {
  suppressPackageStartupMessages(library(CellChat))
}
if (!exists("cco.list")) {
  cco_path <- Sys.getenv(
    "RLN_CELLCHAT_PROCESSED_LIST",
    unset = file.path(datadir, "CellChat", "cco.list_processed.rds")
  )
  if (!file.exists(cco_path)) {
    stop("Provide cco.list or set RLN_CELLCHAT_PROCESSED_LIST for Figure4H.")
  }
  cco.list <- readRDS(cco_path)
}
if (length(cco.list) < 5) {
  stop("Figure4H requires the fifth CellChat object: cco.list[[5]].")
}

p <- netVisual_heatmap(cco.list[[5]], color.heatmap = "Reds")
grDevices::pdf(file.path(figdir, "Figure4H.pdf"), width = 6, height = 5)
ComplexHeatmap::draw(p)
grDevices::dev.off()

####Figure4I####
# Figure 4I: Mac_3-to-myofibroblast rankNet comparison
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Figure4I requires the CellChat package.")
}
if (!"package:CellChat" %in% search()) {
  suppressPackageStartupMessages(library(CellChat))
}
if (!exists("cellchat")) {
  cellchat_path <- Sys.getenv(
    "RLN_CELLCHAT_MERGED_RDS",
    unset = file.path(datadir, "Figure4", "cellchat_merged.rds")
  )
  if (file.exists(cellchat_path)) {
    cellchat <- readRDS(cellchat_path)
  } else {
    cco_path <- Sys.getenv(
      "RLN_CELLCHAT_PROCESSED_LIST",
      unset = file.path(datadir, "CellChat", "cco.list_processed.rds")
    )
    if (!file.exists(cco_path)) {
      stop(
        "Provide merged cellchat, set RLN_CELLCHAT_MERGED_RDS, or set RLN_CELLCHAT_PROCESSED_LIST."
      )
    }
    cco.list <- readRDS(cco_path)
    cellchat <- mergeCellChat(cco.list, add.names = names(cco.list))
  }
}

ranknet_helper <- Sys.getenv("RLN_RANKNET_HELPER", unset = "rankNet_ordered.R")
if (!file.exists(ranknet_helper)) {
  stop("Missing rankNet helper: ", ranknet_helper)
}
source(ranknet_helper)

day <- c("Control", "D3", "D7", "D14", "D28")
day_col <- c("#DC4E32", "#45B6C4", "#259C87", "#3A548C", "#EA987B")
names(day_col) <- day
if (length(methods::slot(cellchat, "netP")) < 5) {
  stop("Figure4I requires at least five CellChat datasets for comparison = 2:5.")
}

p <- rankNet_ordered(
  cellchat,
  mode = "comparison",
  measure = "weight",
  sources.use = "Mac_3",
  color.use = day_col,
  targets.use = "Myofibroblast",
  do.flip = TRUE,
  stacked = TRUE,
  do.stat = TRUE,
  comparison = 2:5,
  return.data = FALSE
)
save_panel("Figure4I.pdf", p, width = 5, height = 4)

####Figure4J####
# Figure 4J: TGFb CellChat bubble plot
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Figure4J requires the CellChat package.")
}
if (!"package:CellChat" %in% search()) {
  suppressPackageStartupMessages(library(CellChat))
}
if (!exists("cellchat")) {
  cco_path <- Sys.getenv(
    "RLN_CELLCHAT_PROCESSED_LIST",
    unset = file.path(datadir, "CellChat", "cco.list_processed.rds")
  )
  if (!file.exists(cco_path)) {
    stop("Provide merged cellchat or set RLN_CELLCHAT_PROCESSED_LIST for Figure4J.")
  }
  cco.list <- readRDS(cco_path)
  cellchat <- mergeCellChat(cco.list, add.names = names(cco.list))
}
if (length(methods::slot(cellchat, "netP")) < 5) {
  stop("Figure4J requires at least five CellChat datasets for comparison = 1:5.")
}

signal_show <- "SPP1"
p <- netVisual_bubble(
  cellchat,
  sources.use = 1:12,
  targets.use = 7:12,
  signaling = "TGFb",
  angle.x = 45,
  remove.isolate = FALSE,
  line.on = TRUE,
  comparison = 1:5,
  return.data = FALSE,
  color.text = c("#DC4E32", "#45B6C4", "#259C87", "#3A548C", "#EA987B")
)

bubble_palette <- colorRampPalette(rev(c(
  "#A40545", "#F46F44", "#FDD985", "#7FCBA4", "#4B65AF"
)))(99)
p <- p +
  scale_colour_gradientn(
    colours = bubble_palette,
    na.value = "white",
    limits = c(quantile(p$data$prob, 0.5, na.rm = TRUE), quantile(p$data$prob, 1, na.rm = TRUE)),
    breaks = c(quantile(p$data$prob, 0.5, na.rm = TRUE), quantile(p$data$prob, 1, na.rm = TRUE)),
    labels = c("min", "max")
  ) +
  facet_wrap(~pathway_name)
save_panel("Figure4J.pdf", p, width = 18, height = 6.5)

####FigS5A####
# Fig S5A: fibroblast monocle3 UMAP
if (exists("Fib")) {
  data <- GetAssayData(Fib, assay = "RNA", slot = "counts")
  cell_metadata <- Fib@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)

  cds <- monocle3::new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )
  cds <- monocle3::preprocess_cds(cds, num_dim = 50)
  cds <- monocle3::reduce_dimension(cds, preprocess_method = "PCA")
  monocle3::plot_cells(cds, reduction_method = "UMAP", color_cells_by = "Fib") +
    ggtitle("cds.umap")

  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(Fib, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed), , drop = FALSE]
  cds@int_colData$reducedDims$UMAP <- int.embed

  col3 <- fib_palette[seq_len(nlevels(Fib))]
  names(col3) <- levels(Fib)
  p <- monocle3::plot_cells(
    cds,
    reduction_method = "UMAP",
    color_cells_by = "Fib",
    show_trajectory_graph = FALSE,
    group_label_size = 8,
    cell_size = 0.5
  ) +
    ggtitle("int.umap") +
    scale_color_manual(values = col3) +
    coord_fixed(ratio = 1.35)
  save_panel("FigS5A.pdf", p, width = 5, height = 5)
}

####FigS5B####
# Fig S5B: fibroblast monocle3 pseudotime UMAP
if (exists("cds")) {
  p <- monocle3::plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    show_trajectory_graph = FALSE,
    cell_stroke = 0.2,
    label_branch_points = FALSE,
    cell_size = 0.8,
    rasterize = FALSE
  ) +
    coord_fixed(ratio = 1.35)
  save_panel("FigS5B.pdf", p, width = 5, height = 5, dpi = 300)
}

####FigS5C####
# Fig S5C: SPP1 CellChat heatmap
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("FigS5C requires the CellChat package.")
}
if (!"package:CellChat" %in% search()) {
  suppressPackageStartupMessages(library(CellChat))
}
if (!exists("cco.list")) {
  cco_path <- Sys.getenv(
    "RLN_CELLCHAT_PROCESSED_LIST",
    unset = file.path(datadir, "CellChat", "cco.list_processed.rds")
  )
  if (!file.exists(cco_path)) {
    stop("Provide cco.list or set RLN_CELLCHAT_PROCESSED_LIST for FigS5C.")
  }
  cco.list <- readRDS(cco_path)
}

object.list <- cco.list
i <- as.integer(Sys.getenv("RLN_S5C_CELLCHAT_INDEX", unset = "5"))
if (is.na(i) || i < 1 || i > length(object.list)) {
  stop("RLN_S5C_CELLCHAT_INDEX must select an existing CellChat object.")
}
pathways.show <- "SPP1"
p <- netVisual_heatmap(
  object.list[[i]],
  signaling = "SPP1",
  color.heatmap = "Reds",
  title.name = paste(pathways.show, "signaling", names(object.list)[i])
)
grDevices::pdf(file.path(figdir, "FigS5C.pdf"), width = 7, height = 5)
ComplexHeatmap::draw(p)
grDevices::dev.off()

####FigS5D####
# Fig S5D: ECM/SPP1 panel
plot_path <- file.path(datadir, "Figure4", "FigS5D.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("FigS5D.pdf", p, width = 7, height = 5)
}

####FigS5E####
# Fig S5E: rankNet panel
plot_path <- file.path(datadir, "Figure4", "FigS5E.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("FigS5E.pdf", p, width = 7, height = 5)
}
