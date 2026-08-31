suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(clusterProfiler)
  library(grid)
  library(nichenetr)
})

# ------------- Figure 3 --------------
# Plotting code for macrophage and T-cell panels.

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
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))

  data.plot$avg.exp.scaled <- avg.exp.scaled
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
      facet_grid(features_group ~ ., scales = "free_y", space = "free_y", switch = "x") +
      theme(panel.spacing = unit(panel.spacing_distance, "lines"), strip.background = element_blank())
  }

  if (return_data) {
    return(list(plot = plot, data = data.plot))
  }
  plot
}

get_sus_scrofa_gobp_list <- function() {
  GO_df_all <- tryCatch(
    msigdbr::msigdbr(species = "Sus scrofa", collection = "C5"),
    error = function(e) msigdbr::msigdbr(species = "Sus scrofa", category = "C5")
  )
  go_subcollection <- if ("gs_subcollection" %in% colnames(GO_df_all)) "gs_subcollection" else "gs_subcat"
  GO_df <- GO_df_all[, c("gs_name", "gene_symbol", "gs_exact_source", go_subcollection)]
  GO_df <- GO_df[GO_df[[go_subcollection]] == "GO:BP", ]
  GO_df$name_id <- paste(tolower(gsub("^GOBP_", "", GO_df$gs_name)), GO_df$gs_exact_source, sep = "_")
  split(GO_df$gene_symbol, GO_df$name_id)
}

mac_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#cab2d6", "#fda111", "#80b1d3", "#8dd3c7")
t_palette <- c("#E47BC0", "#8C87B3", "#E3C4DA", "#89C4D2", "#91D1C2FF", "#89E1F7", "#DFDDF0", "#F27F96", "#B392BD", "#D7F1F6", "#F8C9D2")

mac_path <- file.path(datadir, "Figure3", "Macrophage.rds")
if (file.exists(mac_path)) {
  Macrophage <- readRDS(mac_path)
  mac_levels <- paste0("Mac_", seq_len(6))
  if ("Mac" %in% colnames(Macrophage@meta.data) && all(mac_levels %in% unique(as.character(Macrophage$Mac)))) {
    Macrophage$Mac <- factor(as.character(Macrophage$Mac), levels = mac_levels)
  }
  Idents(Macrophage) <- "Mac"
  mac_palette_use <- rep(mac_palette, length.out = nlevels(Macrophage))
  names(mac_palette_use) <- levels(Macrophage)

  ####Figure3A####
  # Figure 3A: macrophage UMAP
  p <- scplotter::CellDimPlot(
    Macrophage,
    group_by = "Mac",
    reduction = "umap",
    label = TRUE,
    label_insitu = TRUE,
    theme = "theme_blank",
    palcolor = mac_palette_use,
    label_size = 6,
    label_fg = "black",
    shuffle = TRUE,
    label_bg = "black",
    label_bg_r = 0,
    pt_size = 0.1,
    label_repel = TRUE,
    label_repulsion = 80
  )
  save_panel("Figure3A.pdf", p, width = 5, height = 4)

  ####Figure3B####
  # Figure 3B: macrophage marker dotplot
  genes <- list(
    Mac_1 = c("CCL14", "C1QA", "C1QB", "SELENOP", "LYVE1"),
    Mac_2 = c("CCR2", "FCN1", "MEFV", "SELL"),
    Mac_3 = c("MX2", "MAT2A", "SSH2", "DHX34", "ND2"),
    Mac_4 = c("FABP4", "LPL", "FABP5", "MMP19", "CD36"),
    Mac_5 = c("EGLN3", "FLT1", "SERPINE1", "PCOLCE2", "ARG1", "FN1"),
    Cycling = c("CCNA2", "PBK", "KIF2C", "TOP2A", "MKI67")
  )
  p <- yy_Dotplot(
    seuratObj = Macrophage,
    genes = genes,
    group.by = "Mac",
    coord_flip = TRUE,
    gene_pct_cutoff = 0
  )
  p <- p +
    theme(
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank(),
      text = element_text(size = 0),
      panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.2),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.key.size = unit(0.15, "inch"),
      plot.margin = unit(c(0, 0, 0, 0), "char")
    ) +
    labs(x = "", y = "") +
    scale_radius(limits = c(0, 100), range = c(0, 2.5)) +
    coord_flip()
  save_panel("Figure3B.pdf", p, width = 4, height = 7, dpi = 800)

  ####Figure3C####
  # Figure 3C: macrophage composition across time
  p <- SCP::CellStatPlot(
    Macrophage,
    stat.by = "Mac",
    group.by = "orig.ident",
    label = TRUE,
    plot_type = "dot",
    palcolor = mac_palette_use
  ) +
    ylab("") +
    xlab("") +
    theme(legend.position = "none")
  save_panel("Figure3C.pdf", p, width = 3, height = 5)

  ####Figure3D####
  # Figure 3D: pro-/anti-inflammatory score split violin
  pro_inflam_genes <- c("Ccl5", "Ccr7", "Cd40", "Cd86", "Cxcl9", "Cxcl10", "Cxcl11", "Ido1", "Il1a", "Il1b", "Il6", "Irf1", "Irf5", "Kynu")
  anti_inflam_genes <- c("Ccl4", "Ccl13", "Ccl18", "Ccl20", "Ccl22", "Cd276", "Clec7a", "Ctsa", "Ctsb", "Ctsc", "Ctsd", "Fn1", "Il4r", "Irf4", "Lyve1", "Mmp9", "Mmp14", "Mmp19", "Msr1", "Tgfb1", "Tgfb2", "Tgfb3", "Tnfsf8", "Tnfsf12", "Vegfa", "Vegfb", "Vegfc")

  Macrophage <- AddModuleScore(Macrophage, features = list(toupper(pro_inflam_genes)), name = "ProInflamScore")
  Macrophage <- AddModuleScore(Macrophage, features = list(toupper(anti_inflam_genes)), name = "AntiInflamScore")

  meta <- Macrophage@meta.data %>%
    dplyr::select(ProInflamScore1, AntiInflamScore1, Group = Mac) %>%
    mutate(Group = as.factor(Group)) %>%
    pivot_longer(
      cols = c(ProInflamScore1, AntiInflamScore1),
      names_to = "ScoreType",
      values_to = "Score"
    ) %>%
    mutate(ScoreType = ifelse(ScoreType == "ProInflamScore1", "Pro-inflammatory", "Anti-inflammatory"))
  meta$ScoreType <- factor(meta$ScoreType, levels = c("Pro-inflammatory", "Anti-inflammatory"))

  meta_pro <- meta %>% filter(ScoreType == "Pro-inflammatory")
  meta_anti <- meta %>% filter(ScoreType == "Anti-inflammatory")

  custom_colors <- c(
    "Pro-inflammatory" = scales::alpha("#E83828", 0.5),
    "Anti-inflammatory" = scales::alpha("#005A94", 0.5)
  )

  p <- ggplot() +
    gghalves::geom_half_violin(
      data = meta_pro,
      aes(x = Group, y = Score, fill = ScoreType),
      side = "l",
      trim = TRUE,
      scale = "width",
      width = 0.8
    ) +
    gghalves::geom_half_violin(
      data = meta_anti,
      aes(x = Group, y = Score, fill = ScoreType),
      side = "r",
      trim = TRUE,
      scale = "width",
      width = 0.8
    ) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 20),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      legend.margin = margin(t = -5, b = 0),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(size = 6),
        keywidth = 0.8,
        keyheight = 0.6,
        title.position = "top",
        label.position = "right"
      )
  )
  save_panel("Figure3D.pdf", p, width = 8, height = 4.5, dpi = 300)

  mac1_d3_marker_path <- Sys.getenv(
    "RLN_MAC1_D3_MARKERS",
    unset = file.path(datadir, "Figure3", "D3_Mac1_markers.csv")
  )
  if (!file.exists(mac1_d3_marker_path)) {
    mac1_d3_marker_path <- "E:/10X_RLN/new_qc/sub/Macrophage/markers/D3_Mac1_markers.csv"
  }
  if (!file.exists(mac1_d3_marker_path)) {
    stop("Missing marker table for FigS3A/FigS3B/Figure3E: ", mac1_d3_marker_path)
  }
  markers_mac1_d3 <- read.csv(mac1_d3_marker_path, check.names = FALSE)
  markers_mac1_d3 <- markers_mac1_d3[, nzchar(colnames(markers_mac1_d3)), drop = FALSE]
  if (!"gene" %in% colnames(markers_mac1_d3)) {
    markers_mac1_d3$gene <- rownames(markers_mac1_d3)
  }

  need_DEG <- markers_mac1_d3 %>%
    filter(p_val_adj < 0.05) %>%
    transmute(
      log2FoldChange = avg_log2FC,
      pvalue = p_val_adj,
      SYMBOL = gene
    )

  geneList <- need_DEG$log2FoldChange
  names(geneList) <- need_DEG$SYMBOL
  geneList <- sort(geneList, decreasing = TRUE)

  go_list <- get_sus_scrofa_gobp_list()

  ####Figure3E####
  # Figure 3E: GSEA of macrophage pathways
  go_terms <- c(
    "Phagocytosis" = "phagocytosis_GO:0006909",
    "Antigen Processing And Presentation" = "antigen_processing_and_presentation_GO:0019882"
  )
  go_used <- lapply(go_terms, function(term_id) unique(go_list[[term_id]]))
  missing_terms <- names(go_used)[vapply(go_used, is.null, logical(1))]
  if (length(missing_terms) > 0) {
    stop("Missing GO terms for Figure3E: ", paste(missing_terms, collapse = ", "))
  }
  gsea_plots <- lapply(names(go_used), function(term_name) {
    term2gene_single <- data.frame(
      pathway = term_name,
      gene = go_used[[term_name]],
      stringsAsFactors = FALSE
    )
    gsea_cp <- GSEA(
      geneList,
      TERM2GENE = term2gene_single,
      minGSSize = 1,
      pvalueCutoff = 1,
      nPermSimple = 10000,
      verbose = FALSE
    )
    GseaVis::gseaNb(
      gsea_cp,
      geneSetID = term_name,
      subPlot = 2,
      addPval = TRUE,
      pvalX = 0.85,
      pvalY = 0.75,
      nesDigit = 2,
      pDigit = 2,
      pvalSize = 5
    )
  })
  grDevices::pdf(file.path(figdir, "Figure3E.pdf"), width = 5.5, height = 4)
  for (p in gsea_plots) {
    print(p)
  }
  grDevices::dev.off()

  ####Figure3G####
  # Figure 3G: macrophage functional pathway scores
  go_description <- c(
    "phagocytosis_engulfment_GO:0006911",
    "antigen_processing_and_presentation_of_endogenous_antigen_GO:0019883",
    "positive_regulation_of_alpha_beta_t_cell_activation_GO:0046635",
    "cholesterol_storage_GO:0010878",
    "lipid_transport_across_blood_brain_barrier_GO:1990379"
  )
  df <- Macrophage@meta.data[, c("orig.ident", "Mac"), drop = FALSE]
  for (item in go_description) {
    descrip <- sub("_GO:\\d+$", "", item)
    df[, descrip] <- AddModuleScore(
      Macrophage,
      go_list[item],
      name = "addscore",
      search = FALSE,
      assay = "RNA"
    )$addscore1
  }
  df_long <- df %>%
    pivot_longer(
      cols = -c(orig.ident, Mac),
      names_to = "sig",
      values_to = "score"
    )
  df_long$Mac <- factor(df_long$Mac, levels = levels(Macrophage))

  p <- ggplot(df_long, aes(
    x = Mac,
    y = score,
    color = Mac,
    fill = Mac
  )) +
    geom_boxplot(outlier.color = NA, lwd = 0.3, notch = TRUE) +
    scale_color_manual(values = mac_palette_use) +
    scale_fill_manual(values = scales::alpha(mac_palette_use, 0.2)) +
    labs(x = "", y = "") +
    facet_wrap(~sig, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      text = element_text(size = 8),
      plot.title = element_text(size = 7, hjust = 0, face = "plain"),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "char"),
      strip.background = element_blank()
    )
  save_panel("Figure3G.pdf", p, width = 1.75, height = 5)

}

t_path <- file.path(datadir, "Figure3", "T_cell.rds")
if (file.exists(t_path)) {
  T_cell <- readRDS(t_path)
  t_group_col <- if ("T_type" %in% colnames(T_cell@meta.data)) "T_type" else "T"
  Idents(T_cell) <- t_group_col
  t_palette_use <- rep(t_palette, length.out = nlevels(T_cell))
  names(t_palette_use) <- levels(T_cell)

  ####Figure3H####
  # Figure 3H: T-cell UMAP
  p <- scplotter::CellDimPlot(
    T_cell,
    group_by = t_group_col,
    reduction = "umap",
    label = TRUE,
    label_insitu = TRUE,
    theme = "theme_blank",
    palcolor = t_palette_use,
    label_size = 6,
    label_fg = "black",
    shuffle = TRUE,
    label_bg = "black",
    label_bg_r = 0,
    pt_size = 0.1,
    label_repel = TRUE,
    label_repulsion = 80
  )
  save_panel("Figure3H.pdf", p, width = 5, height = 4)

  ####Figure3I####
  # Figure 3I: T subtype GO enrichment barplot
  t_enrich_path <- file.path(datadir, "Figure3", "Tcell_GO_terms.csv")
  if (file.exists(t_enrich_path)) {
    select_result <- read.csv(t_enrich_path, check.names = FALSE)
    if ("select" %in% colnames(select_result)) {
      select_result <- select_result %>% filter(select == 1)
    }
    fig3i_clusters <- levels(T_cell)[c(1, 2, 4)]
    fig3i_clusters <- fig3i_clusters[!is.na(fig3i_clusters)]
    if (all(fig3i_clusters %in% unique(as.character(select_result$cluster)))) {
      select_result <- select_result %>% filter(cluster %in% fig3i_clusters)
      select_result$cluster <- factor(select_result$cluster, levels = fig3i_clusters)
    } else {
      select_result$cluster <- factor(select_result$cluster, levels = unique(select_result$cluster))
    }
    select_result$Description <- factor(
      select_result$Description,
      levels = select_result$Description[order(select_result$qvalue)] %>% unique()
    )
    t_enrich_palette <- t_palette_use[levels(select_result$cluster)]
    if (any(is.na(t_enrich_palette))) {
      t_enrich_palette <- rep(t_palette, length.out = nlevels(select_result$cluster))
      names(t_enrich_palette) <- levels(select_result$cluster)
    }

    p <- ggplot(select_result, aes(x = Description, y = -log(qvalue), fill = cluster)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = t_enrich_palette) +
      facet_grid(rows = vars(cluster), scales = "free", space = "free") +
      theme_bw() +
      theme(
        panel.grid = element_blank()
      ) +
      xlab("") +
      ylab("-log10(FDR)") +
      theme(
        strip.text.x = element_text(colour = "black", angle = 360, size = 12, hjust = 0, face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.line = element_line(size = 1, colour = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      coord_flip()
    save_panel("Figure3I.pdf", p, width = 10, height = 8)
  }
}

####Figure3J####
# Figure 3J: NicheNet signaling network
nichenet_dir <- Sys.getenv("RLN_NICHENET_DIR", unset = file.path(datadir, "Figure3", "NicheNet"))
if (!dir.exists(nichenet_dir)) {
  nichenet_dir <- "E:/10X_RLN/new_qc/Nichenet"
}
nichenet_files <- file.path(
  nichenet_dir,
  c(
    "lr_network_human_21122021.rds",
    "signaling_network_human_21122021.rds",
    "weighted_networks_nsga2r_final.rds",
    "gr_network_human_21122021.rds",
    "ligand_tf_matrix_nsga2r_final.rds"
  )
)
if (all(file.exists(nichenet_files))) {
  lr_network <- readRDS(file.path(nichenet_dir, "lr_network_human_21122021.rds"))
  sig_network <- readRDS(file.path(nichenet_dir, "signaling_network_human_21122021.rds"))
  weighted_networks <- readRDS(file.path(nichenet_dir, "weighted_networks_nsga2r_final.rds"))
  gr_network <- readRDS(file.path(nichenet_dir, "gr_network_human_21122021.rds"))
  ligand_tf_matrix <- readRDS(file.path(nichenet_dir, "ligand_tf_matrix_nsga2r_final.rds"))

  ligands_oi <- c("IL17C", "IL1A", "IFNA17", "IFNB1", "SPARC")
  targets_oi <- c("SPP1", "SPI1")

  active_signaling_network <- get_ligand_signaling_path(
    ligands_all = ligands_oi,
    targets_all = targets_oi,
    weighted_networks = weighted_networks,
    ligand_tf_matrix = ligand_tf_matrix,
    top_n_regulators = 4,
    minmax_scaling = TRUE
  )

  graph_min_max <- diagrammer_format_signaling_graph(
    signaling_graph_list = active_signaling_network,
    ligands_all = ligands_oi,
    targets_all = targets_oi,
    sig_color = "indianred",
    gr_color = "steelblue"
  )

  graph_svg <- DiagrammeRsvg::export_svg(
    DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph")
  )
  figure3j_svg <- file.path(figdir, "Figure3J.svg")
  writeLines(graph_svg, figure3j_svg)
  rsvg::rsvg_pdf(figure3j_svg, file.path(figdir, "Figure3J.pdf"))
}

if (exists("T_cell")) {
  ####Figure3K####
  # Figure 3K: type II IFN production score
  go_list <- get_sus_scrofa_gobp_list()
  path.name <- "type_ii_interferon_production"
  target <- list()
  target[[path.name]] <- go_list[["type_ii_interferon_production_GO:0032609"]]
  if (!is.null(target[[path.name]])) {
    df <- T_cell@meta.data[, c("orig.ident", t_group_col), drop = FALSE]
    df$Group <- factor(df[[t_group_col]], levels = levels(T_cell))
    kk <- if (exists("AddModuleScore_UCell", mode = "function")) {
      AddModuleScore_UCell(T_cell, target, name = "addscore")@meta.data
    } else {
      UCell::AddModuleScore_UCell(T_cell, target, name = "addscore")@meta.data
    }
    df$addscore <- kk[, ncol(kk)]

    sig_mean <- df %>%
      dplyr::summarise(mean_score = mean(addscore, na.rm = TRUE))
    stat_df <- df %>%
      group_by(Group) %>%
      dplyr::summarise(
        Q1 = quantile(addscore, 0.25, na.rm = TRUE),
        Q3 = quantile(addscore, 0.75, na.rm = TRUE),
        median = median(addscore, na.rm = TRUE),
        .groups = "drop"
      )
    stat_df$Group <- factor(stat_df$Group, levels = levels(df$Group))

    comparison_idx <- list(c(2, 4), c(3, 4))
    my_comparisons <- lapply(comparison_idx, function(idx) levels(df$Group)[idx])
    my_comparisons <- Filter(function(x) length(x) == 2 && all(!is.na(x)), my_comparisons)

    p <- ggplot(df, aes(x = Group, y = addscore, fill = Group, color = Group)) +
      geom_violin(
        scale = "width",
        adjust = 1,
        trim = TRUE
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
        aes(x = Group, ymin = Q1, ymax = Q3),
        width = 0.1,
        color = "black",
        inherit.aes = FALSE,
        linewidth = 0.8
      ) +
      geom_segment(
        data = stat_df,
        inherit.aes = FALSE,
        aes(
          x = as.numeric(Group) - 0.05,
          xend = as.numeric(Group) + 0.05,
          y = median,
          yend = median
        ),
        color = "black",
        linewidth = 0.8
      ) +
      labs(title = path.name) +
      theme_classic() +
      ylab("Gene score") +
      xlab("") +
      scale_color_manual(values = t_palette_use) +
      scale_fill_manual(values = t_palette_use)
    if (length(my_comparisons) > 0) {
      p <- p +
        ggpubr::stat_compare_means(
          comparisons = my_comparisons,
          method = "wilcox.test",
          label = "p.signif"
        )
    }
    p <- p +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 22),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 18, hjust = 0.5, vjust = 1, color = "black")
      )
    save_panel("Figure3K.pdf", p, width = 8, height = 6, dpi = 300)
  }
}

if (exists("Macrophage")) {
  ####FigS3A####
  # Fig S3A: Mac_1 D3 vs Control volcano plot
  p <- plotthis::VolcanoPlot(
    markers_mac1_d3,
    x = "avg_log2FC",
    y = "p_val_adj",
    y_cutoff_name = "none",
    label_by = "gene",
    labels = c(
      "LYVE1", "COLEC12", "C1QC", "IL1A", "IL1RAP",
      "SELENOP", "IL1B", "TNFAIP6", "LGALS3", "MARCO"
    ),
    label_size = 7,
    pt_size = 0.5
  )
  save_panel("FigS3A.pdf", p, width = 5, height = 5, dpi = 300)

  ####FigS3B####
  # Fig S3B: Type II interferon mediated signaling pathway GSEA
  go_used_s3b <- list(
    "Type II Interferon Mediated Signaling Pathway" =
      unique(go_list[["type_ii_interferon_mediated_signaling_pathway_GO:0060333"]])
  )
  if (is.null(go_used_s3b[[1]])) {
    stop("Missing GO term for FigS3B: type_ii_interferon_mediated_signaling_pathway_GO:0060333")
  }
  go_used_s3b_df <- stack(go_used_s3b)
  colnames(go_used_s3b_df) <- c("gene", "term_id")

  gsea_s3b <- GSEA(
    geneList,
    TERM2GENE = go_used_s3b_df[, c("term_id", "gene")],
    minGSSize = 1,
    pvalueCutoff = 1,
    nPermSimple = 10000,
    verbose = FALSE
  )
  p <- GseaVis::gseaNb(
    gsea_s3b,
    geneSetID = names(go_used_s3b),
    subPlot = 2,
    addPval = TRUE,
    pvalX = 0.85,
    pvalY = 0.75,
    nesDigit = 2,
    pDigit = 2,
    pvalSize = 5,
    ncol = 1
  )
  save_panel("FigS3B.pdf", p, width = 5.5, height = 4)

  ####FigS3C####
  # Fig S3C: PROGENy pathway activity plot
  options(RCurlOptions = list(ftp.use.epsv = FALSE))
  progeny_net_path <- Sys.getenv(
    "RLN_PROGENY_NET",
    unset = file.path(datadir, "Figure3", "progeny_net.rds")
  )
  if (!file.exists(progeny_net_path)) {
    progeny_net_path <- "E:/10X_RLN/new_qc/sub/T_cell/progeny_net.rds"
  }
  if (!file.exists(progeny_net_path)) {
    stop("Missing PROGENy network for FigS3C: ", progeny_net_path)
  }

  net <- readRDS(progeny_net_path)
  mat <- as.matrix(Macrophage[["Homo_RNA"]]$data)
  acts <- decoupleR::run_mlm(
    mat = mat,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = 5
  )
  Macrophage[["pathwaysmlm"]] <- acts %>%
    pivot_wider(
      id_cols = "source",
      names_from = "condition",
      values_from = "score"
    ) %>%
    tibble::column_to_rownames("source") %>%
    Seurat::CreateAssayObject(.)

  DefaultAssay(Macrophage) <- "pathwaysmlm"
  Macrophage <- ScaleData(Macrophage)
  Macrophage[["pathwaysmlm"]]$data <- Macrophage[["pathwaysmlm"]]$scale.data

  progeny_scores_df <- as.data.frame(t(GetAssayData(
    Macrophage,
    slot = "scale.data",
    assay = "pathwaysmlm"
  ))) %>%
    tibble::rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell)
  CellsClusters <- data.frame(
    Cell = colnames(Macrophage),
    CellType = as.character(Macrophage@meta.data$Mac),
    stringsAsFactors = FALSE
  )
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

  summarized_progeny_scores <- progeny_scores_df %>%
    dplyr::group_by(Pathway, CellType) %>%
    dplyr::summarise(avg = mean(Activity), std = sd(Activity), .groups = "drop")
  sub_progeny_scores <- summarized_progeny_scores %>%
    filter(Pathway %in% c("NFkB", "TNFa", "TGFb", "Hypoxia", "VEGF"))

  cluster.lineages <- paste0("Mac_", c(1:5))
  sub_progeny_scores <- sub_progeny_scores %>%
    filter(CellType %in% cluster.lineages)
  sub_progeny_scores$CellType <- factor(sub_progeny_scores$CellType, levels = rev(cluster.lineages))
  sub_progeny_scores$Pathway <- factor(sub_progeny_scores$Pathway, levels = c("NFkB", "TNFa", "TGFb", "Hypoxia", "VEGF"))
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
  progeny_palette <- c("#ADD8E6", "#1E78B4", "#B2E281", "#2CA02C", "#FFB347")

  p <- ggplot(
    sub_progeny_scores,
    aes(x = CellType, y = avg, colour = group)
  ) +
    geom_segment(
      data = line_df,
      aes(
        x = CellType,
        xend = CellType,
        y = y1,
        yend = y2,
        colour = group
      ),
      linetype = "solid",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(size = 2.4) +
    geom_hline(
      data = groupMean,
      aes(yintercept = MeanV),
      linetype = "dashed"
    ) +
    coord_flip() +
    scale_colour_manual(values = progeny_palette) +
    facet_grid(group ~ Pathway, scales = "free") +
    theme_gray(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7),
      strip.background = element_blank()
    )
  save_panel("FigS3C.pdf", p, width = 6, height = 3.5)

  ####FigS3D####
  # Fig S3D: macrophage GO enrichment dotplot
  mac_go_path <- Sys.getenv(
    "RLN_MAC_GO_ANNO",
    unset = file.path(datadir, "Figure3", "fig3_Hs_DEG_GO_anno.csv")
  )
  if (!file.exists(mac_go_path)) {
    mac_go_path <- "E:/10X_RLN/new_qc/sub/Macrophage/markers/fig3_Hs_DEG_GO_anno.csv"
  }
  if (!file.exists(mac_go_path)) {
    stop("Missing GO enrichment table for FigS3D: ", mac_go_path)
  }
  filtered_go <- read.csv(mac_go_path, check.names = FALSE) %>%
    filter(select == 1) %>%
    mutate(
      Description = paste0(
        toupper(substr(Description, 1, 1)),
        substr(Description, 2, nchar(Description))
      )
    )
  filtered_go$Cluster <- factor(filtered_go$Cluster, levels = levels(Macrophage))

  p <- scplotter::EnrichmentPlot(
    filtered_go,
    plot_type = "comparison",
    group_by = "Cluster",
    top_term = 10,
    aspect.ratio = 2,
    x_text_angle = 45
  ) +
    scale_fill_gradientn(
      colours = brewer.pal(11, "Spectral") %>% rev(),
      limits = c(0, 5),
      oob = scales::squish,
      name = "-log10(p.adjust)"
    ) +
    scale_y_discrete(limits = rev)
  save_panel("FigS3D.pdf", p, width = 9, height = 10, dpi = 300)

  ####FigS3E####
  # Fig S3E: D28 macrophage SCENIC regulon specificity score
  scenic_loom_path <- Sys.getenv(
    "RLN_MAC_D28_SCENIC_LOOM",
    unset = file.path(datadir, "Figure3", "input_counts_D28_scenic.loom")
  )
  if (!file.exists(scenic_loom_path)) {
    scenic_loom_path <- "E:/10X_RLN/new_qc/sub/Macrophage/scenic/input_counts_D28_scenic.loom"
  }
  if (!file.exists(scenic_loom_path)) {
    stop("Missing SCENIC loom file for FigS3E: ", scenic_loom_path)
  }
  if (!requireNamespace("SCENIC", quietly = TRUE)) {
    stop("Package SCENIC is required for FigS3E.")
  }
  if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
    stop("Package SCopeLoomR is required for FigS3E.")
  }
  suppressPackageStartupMessages({
    library(SCENIC)
    library(SCopeLoomR)
  })

  loom <- open_loom(scenic_loom_path)
  exprMat <- get_dgem(loom)
  exprMat_log <- log2(exprMat + 1)
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  regulonAucThresholds <- get_regulon_thresholds(loom)
  close_loom(loom)

  sce <- subset(Macrophage, orig.ident == "D28")
  sub_regulonAUC <- regulonAUC[, match(rownames(sce@meta.data), colnames(regulonAUC))]
  if (!identical(colnames(sub_regulonAUC), rownames(sce@meta.data))) {
    stop("FigS3E SCENIC AUC columns do not match D28 macrophage metadata rows.")
  }

  cellTypes <- data.frame(
    row.names = rownames(sce@meta.data),
    celltype = sce$Mac
  )
  rss <- calcRSS(AUC = getAUC(sub_regulonAUC), sce$Mac)
  rss <- na.omit(rss)
  rss <- rss[, paste0("Mac_", c(1:6))]
  p <- plotRSS(
    rss,
    order_rows = TRUE,
    labelsToDiscard = c("Mac_1", "Mac_2", "Mac_3", "Mac_6"),
    zThreshold = 1.1
  )
  save_panel("FigS3E.pdf", p, width = 7, height = 5)
}

####FigS4A####
# Fig S4A: CellChat differential interaction weight
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("FigS4A requires the CellChat package.")
}
if (!"package:CellChat" %in% search()) {
  suppressPackageStartupMessages(library(CellChat))
}
if (!exists("cellchat")) {
  cellchat_path <- Sys.getenv(
    "RLN_CELLCHAT_MERGED_RDS",
    unset = file.path(datadir, "CellChat", "cellchat_merged.rds")
  )
  if (file.exists(cellchat_path)) {
    cellchat <- readRDS(cellchat_path)
  } else {
    cco_path <- Sys.getenv(
      "RLN_CELLCHAT_PROCESSED_LIST",
      unset = file.path(datadir, "CellChat", "cco.list_processed.rds")
    )
    if (!file.exists(cco_path)) {
      stop("Provide merged cellchat, set RLN_CELLCHAT_MERGED_RDS, or set RLN_CELLCHAT_PROCESSED_LIST.")
    }
    cco.list <- readRDS(cco_path)
    cellchat <- mergeCellChat(cco.list, add.names = names(cco.list))
  }
}
if (length(methods::slot(cellchat, "net")) < 2) {
  stop("FigS4A requires at least two CellChat datasets.")
}

p <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
save_panel("FigS4A.pdf", p, width = 7, height = 5)

####FigS4B####
# Fig S4B: T-cell marker heatmap
if (exists("T_cell")) {
  t_marker_path <- Sys.getenv(
    "RLN_T_MARKERS",
    unset = file.path(datadir, "Figure3", "findallmarkers11.csv")
  )
  if (!file.exists(t_marker_path)) {
    t_marker_path <- "E:/10X_RLN/new_qc/sub/T_cell/markers/findallmarkers11.csv"
  }
  if (!file.exists(t_marker_path)) {
    stop("Missing marker table for FigS4B: ", t_marker_path)
  }
  top10 <- read.csv(t_marker_path) %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC)
  top10$cluster <- factor(top10$cluster, levels = levels(T_cell))
  top10 <- top10[order(top10$cluster), ]

  t_heatmap_levels <- levels(T_cell)[c(1, 2, 4, 3, 5, 6, 7)]
  t_heatmap_levels <- t_heatmap_levels[!is.na(t_heatmap_levels)]
  t_heatmap_colors <- t_palette[c(1, 2, 4, 3, 5, 6, 7)][seq_along(t_heatmap_levels)]
  names(t_heatmap_colors) <- t_heatmap_levels
  T_cell[[t_group_col]] <- factor(T_cell[[t_group_col, drop = TRUE]], levels = t_heatmap_levels)
  Idents(T_cell) <- t_group_col

  scop::show_palettes()
  bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
    cell_type = levels(T_cell),
    col = list(cell_type = t_heatmap_colors),
    annotation_height = grid::unit(6, "mm")
  )
  ht <- scop::GroupHeatmap(
    srt = T_cell,
    label_size = 10,
    show_column_names = TRUE,
    column_names_side = "bottom",
    width = 4,
    height = 6,
    features = c(top10$gene %>% unique(), "IL23R", "IL17A", "LOC102160313"),
    features_label = c(
      "CD40LG", "CCR4", "IL7R", "IL6ST", "FURIN", "CCR7",
      "FOXP3", "IKZF4", "CTLA4", "IL2RA", "CCR8", "TNFRSF18", "IKZF2",
      "LEF1", "TCF7", "CD27", "SELL", "ID3", "S1PR1",
      "GZMK", "GZMB", "CCL4", "FCGR3A", "PRF1",
      "NCR1", "GNLY", "KLRD1", "NKG7",
      "BLK", "IL23R", "SOX13", "IL17A", "CCR6", "LOC102160313"
    ),
    group.by = t_group_col,
    group_palcolor = list(t_heatmap_colors),
    column_names_rot = 30,
    column_title_side = "bottom",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_palcolor = colorRampPalette(c("#1269B2", "#FFFDF2", "#BC4A4A"))(40)
  )
  save_panel("FigS4B.pdf", ht$plot, width = 8, height = 12, dpi = 800)
}

####FigS4C####
# Fig S4C: T-cell G2M panel
if (exists("T_cell")) {
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  T_cell <- CellCycleScoring(
    T_cell,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )
  p <- scplotter::FeatureStatPlot(
    T_cell,
    features = c("G2M.Score"),
    palcolor = t_palette_use,
    ident = t_group_col,
    plot_type = "bar",
    facet_scales = "free_y"
  )
  save_panel("FigS4C.pdf", p, width = 8, height = 5)
}

####FigS4D####
# Fig S4D: T-cell GO enrichment lollipop plot
t_go_path <- Sys.getenv(
  "RLN_T_GO",
  unset = file.path(datadir, "Figure3", "Tcell_GO_terms.csv")
)
if (!file.exists(t_go_path)) {
  t_go_path <- "E:/10X_RLN/new_qc/sub/T_cell/enrich/enrich_go.csv"
}
if (!file.exists(t_go_path)) {
  stop("Missing GO enrichment table for FigS4D: ", t_go_path)
}
select_result <- read.csv(t_go_path) %>%
  filter(select == 1) %>%
  group_by(cluster) %>%
  arrange(cluster, desc(log(qvalue)))
p <- scplotter::EnrichmentPlot(
  select_result,
  plot_type = "lollipop",
  top_term = 10,
  facet_by = "cluster",
  facet_nrow = 3,
  aspect.ratio = 0.4
)
save_panel("FigS4D.pdf", p, width = 10, height = 6)

####FigS4E####
# Fig S4E: PRF1 panel
plot_path <- file.path(datadir, "Figure3", "FigS4E.rds")
if (file.exists(plot_path)) {
  p <- readRDS(plot_path)
  save_panel("FigS4E.pdf", p, width = 7, height = 5)
}
