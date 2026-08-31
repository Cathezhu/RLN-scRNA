suppressPackageStartupMessages({
  library(CellChat)
})

# ------------- CellChat --------------
# Input: a named list of CellChat objects created for each time point.
# The objects must retain their expression data and cell-type identities.

cellchat_datadir <- Sys.getenv("RLN_CELLCHAT_DATA", unset = "./data/CellChat")
cellchat_outdir <- Sys.getenv("RLN_CELLCHAT_OUTPUT", unset = "./figures/CellChat")
dir.create(cellchat_outdir, recursive = TRUE, showWarnings = FALSE)

day <- c("Control", "D3", "D7", "D14", "D28")
day_col <- c("#DC4E32", "#45B6C4", "#259C87", "#3A548C", "#EA987B")
names(day_col) <- day

####CellChatInput####
if (!exists("cco.list")) {
  cco_path <- Sys.getenv(
    "RLN_CELLCHAT_INPUT",
    unset = file.path(cellchat_datadir, "cco.list.rds")
  )
  if (!file.exists(cco_path)) {
    stop(
      "Provide a named CellChat object list as cco.list, or set RLN_CELLCHAT_INPUT to its RDS file."
    )
  }
  cco.list <- readRDS(cco_path)
}
if (!is.list(cco.list) || length(cco.list) == 0) {
  stop("cco.list must be a non-empty named list of CellChat objects.")
}
if (is.null(names(cco.list)) || any(!nzchar(names(cco.list)))) {
  if (length(cco.list) == length(day)) {
    names(cco.list) <- day
  } else {
    stop("Name cco.list with dataset/time-point names before running CellChat.")
  }
}

####CellChatInference####
for (i in seq_along(cco.list)) {
  cco.list[[i]]@idents <- droplevels(cco.list[[i]]@idents)
  print(levels(cco.list[[i]]@idents))
}

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(
  CellChatDB,
  search = c("Secreted Signaling", "Cell-Cell Contact"),
  key = "annotation"
)

for (i in seq_along(cco.list)) {
  future::plan("sequential")
  message("1. New dataset starts: No.", i, "!!")

  cellchat_obj <- cco.list[[i]]
  cellchat_obj@DB <- CellChatDB.use
  groupSize <- as.numeric(table(cellchat_obj@idents))
  cellchat_obj <- subsetData(cellchat_obj)

  message("2. Processing identifyOverExpressedGenes!!!")
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  message("3. Processing identifyOverExpressedInteractions!!!")
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  message("4. Processing projectData!!!")
  cellchat_obj <- projectData(cellchat_obj, PPI.human)
  message("5. Processing computeCommunProb!!!")
  cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean", raw.use = FALSE)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  message("6. Processing computeCommunProbPathway!!!")
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)

  cco.list[[i]] <- cellchat_obj
}

gc()
saveRDS(cco.list, file.path(cellchat_outdir, "cco.list_processed.rds"))

####SPP1Circle####
pathways.show <- "SPP1"
circle_list <- cco.list[seq_len(min(3, length(cco.list)))]
weight.max <- getMaxWeight(
  circle_list,
  slot.name = "netP",
  attribute = pathways.show
)
png(
  file.path(cellchat_outdir, paste0(pathways.show, "_circle.png")),
  height = 4500,
  width = 3000 * length(circle_list),
  res = 300
)
par(mfrow = c(1, length(circle_list)), xpd = TRUE)
for (i in seq_along(circle_list)) {
  netVisual_aggregate(
    circle_list[[i]],
    thresh = 0.0001,
    signaling = pathways.show,
    sources.use = "Mac_3",
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, names(circle_list)[i])
  )
}
dev.off()

####DifferentialInteraction####
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list))
if (length(cco.list) >= 2) {
  comparisons <- lapply(seq_len(length(cco.list) - 1), function(i) c(i, i + 1))
  png(
    file.path(cellchat_outdir, "differential_interaction.png"),
    height = 2000,
    width = 1800 * length(comparisons),
    res = 550
  )
  par(mfrow = c(1, length(comparisons)), xpd = TRUE)
  for (comparison in comparisons) {
    color.use <- day_col[names(cco.list)[comparison]]
    if (any(is.na(color.use))) {
      color.use <- c("#DC4E32", "#45B6C4")
    }
    netVisual_diffInteraction(
      cellchat,
      comparison = comparison,
      targets.use = "Endothelial_cell",
      weight.scale = TRUE,
      measure = "weight",
      remove.isolate = FALSE,
      vertex.label.cex = 0.9,
      margin = 0.1,
      color.use = color.use
    )
  }
  dev.off()
}
