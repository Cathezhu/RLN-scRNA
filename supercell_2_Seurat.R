+supercell_2_Seurat <- function(SC.GE, SC, fields = c(),
                               var.genes = NULL,
                               do.preproc = TRUE,
                               is.log.normalized = TRUE,
                               do.center = TRUE,
                               do.scale = TRUE,
                               N.comp = NULL,
                               output.assay.version = "v4"){
  N.c <- ncol(SC.GE)
  if(is.null(SC$supercell_size)){
    warning(paste0("supercell_size field of SC is missing, size of all super-cells set to 1"))
    supercell_size <- rep(1, N.c)
  } else {
    supercell_size <- SC$supercell_size
  }

  if(length(supercell_size) != N.c){
    stop(paste0("length of SC$supercell_size has to be the same as number of super-cells ", N.c))
  }

  ## Name all cells to create Seurat Object
  if(is.null(colnames(SC.GE))){
    colnames(SC.GE) <- as.character(1:N.c)
  }
  counts <- SC.GE

  ## If fields is numerical, map them to names
  if(is.numeric(fields)){
    fields <- names(SC)[fields]
  }

  ## Keep only available fiedls
  fields <- intersect(fields, names(SC))

  if(length(fields) > 0){
    SC.fields <- SC[fields]
  } else {
    SC.fields <- NULL
  }

  ## Keep only fields that are specific to cells
  SC.field.length <- lapply(SC.fields, length)
  SC.fields       <- SC.fields[which(SC.field.length == N.c)]

  meta     <- data.frame(size = supercell_size, row.names = colnames(SC.GE), stringsAsFactors = FALSE)

  if(length(SC.fields) > 0){
    meta <- cbind(meta, SC.fields)
  }
  m.seurat <- Seurat::CreateSeuratObject(counts = SC.GE, meta.data = meta)


  if(!do.preproc) return(m.seurat)

  ## Data preprocessing (optional, but recommended)
  ## Normalize data, so Seurat does not generate warning
  m.seurat <- Seurat::NormalizeData(m.seurat)
  print("Done: NormalizeData")
  ## If SC.GE is log-normalized gene expression, than field data has to be rewritten
  if(is.log.normalized){
    print("Doing: data to normalized data")
    m.seurat[["RNA"]]$data <- m.seurat[["RNA"]]$counts
  }

  ## Sample-weighted scaling
  if(length(unique(meta$size)) > 1){
    print("Doing: weighted scaling")
    m.seurat[["RNA"]]$scale.data <- t(as.matrix(corpcor::wt.scale(Matrix::t((m.seurat[["RNA"]]$data)),
                                                                    w = meta$size,
                                                                    center = do.center,
                                                                    scale = do.scale)))
    print("Done: weighted scaling")

  } else {
    print("Doing: unweighted scaling")
    m.seurat <- Seurat::ScaleData(m.seurat)
    print("Done: unweighted scaling")
  }

  m.seurat@assays$RNA@misc[["scale.data.weighted"]] <- m.seurat[["RNA"]]$scale.data
+
  if(is.null(var.genes)){
    var.genes <- sort(SC$genes.use)
  }


  if(is.null(N.comp)) N.comp <- min(50, ncol(m.seurat[["RNA"]]$counts)-1)

  Seurat::VariableFeatures(m.seurat) <- var.genes
  m.seurat <- Seurat::RunPCA(m.seurat, verbose = F, npcs = max(N.comp))
  m.seurat@reductions$pca_seurat <- m.seurat@reductions$pca

  my_pca <- supercell_prcomp(X = Matrix::t(SC.GE[var.genes, ]), genes.use = var.genes,
                             fast.pca = TRUE,
                             supercell_size = meta$supercell_size,
                             k = dim(m.seurat@reductions$pca_seurat)[2],
                             do.scale = do.scale, do.center = do.center)

  dimnames(my_pca$x) <- dimnames(m.seurat@reductions$pca_seurat)
  m.seurat@reductions$pca@cell.embeddings  <- my_pca$x
  m.seurat@reductions$pca@feature.loadings <- my_pca$rotation
  m.seurat@reductions$pca@stdev            <- my_pca$sdev

  m.seurat@reductions$pca_weighted         <- m.seurat@reductions$pca

  ## Super-cell network:
  ## 1) create graph field
  m.seurat            <- Seurat::FindNeighbors(m.seurat, compute.SNN = TRUE, verbose = TRUE)

  ## 2) add self-loops to our super-cell graph to indicate super-cell size (does not work, as Seurat removes loops...)

  if(!is.null(SC$graph.supercells)){
    # SC$graph.supercells <- igraph::add_edges(SC$graph.supercells, edges = rep(1:N.c, each = 2), weight = supercell_size)
    adj.mtx             <- igraph::get.adjacency(SC$graph.supercells, attr = "weight")

    ## 3) replace generated Seurat network with the super-cell network
    m.seurat@graphs$RNA_nn@i                <- adj.mtx@i
    m.seurat@graphs$RNA_nn@p                <- adj.mtx@p
    m.seurat@graphs$RNA_nn@Dim              <- adj.mtx@Dim
    m.seurat@graphs$RNA_nn@x                <- adj.mtx@x
    m.seurat@graphs$RNA_nn@factors          <- adj.mtx@factors

    m.seurat@graphs$RNA_super_cells         <- m.seurat@graphs$RNA_nn
  } else {
    warning("Super-cell graph was not found in SC object, no super-cell graph was added to Seurat object")
  }
  if(as.character(packageVersion("Seurat")) >= "5.0.0" & output.assay.version == "v5"){
    m.seurat[["RNA"]] <- as(object = m.seurat[["RNA"]], Class = "Assay5")
  }

  return(m.seurat)
}
# Seurat conversion helper used by the SuperCell workflow.
