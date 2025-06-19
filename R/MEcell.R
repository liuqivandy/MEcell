#'MEcell

#' @description Leverage the microenvironment profile to improve spatial cell clustering
#'
#' @param obj A Seurat object
#' @param assay Assay to use in clustering. Default is 'NULL', that is the Default Assay of obj
#' @param k_spatial The number of spatial neighbors used to compute microenvironment profile. Default is 16
#' @param k_nn The number of nearest neighbors. Default is 20
#' @param usepca whether to use PCA scores to compute microenvironment profile. Default is FALSE. Set to True for whole transcriptome platform for speed.
#' @param K_adp whether to use adaptive strategy. Default is FALSE.
#' @param delta The parameter controlling the sensitivity to local microenvironment change. Default is 0.5,only valid when k_adp=T
#' @param prune.SNN The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15
#'
#' @return A Seurat object with a microenvironment-refined SNN graph (named as "MEcell") stored in the respective slot.
#'
#' library(MEcell)
#' obj<-readRDS(url("https://www.dropbox.com/s/f5khi1zperqkybg/Mouse_Brain_Serial_Section1_SagittalPosterior_rawimage_object.rds?dl=1"))
#' obj <- FindClusters(obj,graph.name="MEcell")

#' @export



MEcell<-function(obj,assay=NULL,k_spatial=16,k_nn=20,usepca=F,K_adp=F,delta=0.5,prune.SNN=1/15){

  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
    message("Using assay: ", assay)
  }

  # Check if assay exists
  if (!(assay %in% names(obj@assays))) {
    stop("Assay not found: ", assay)
  }

  # Get the assay object
  assay_obj <- obj[[assay]]

  # Determine if Seurat v5 (Assay5) or v4 (Assay)
  is_v5 <- inherits(assay_obj, "Assay5")

  # Check and retrieve 'data'
  has_data <- FALSE

  if (is_v5) {
    has_data <- "data" %in% names(assay_obj@layers)
  } else {
    has_data <- !is.null(assay_obj@data) && ncol(assay_obj@data) > 0
  }

  # Create 'data' if missing
  if (!has_data) {
    message("'data' not found. Running NormalizeData()...")
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
    obj<-ScaleData(obj,features=rownames(obj))
  }


  if (!"pca" %in% names(obj@reductions)){

    message("PCA not found — running PCA now...")
    obj <- RunPCA(obj, verbose = FALSE)
  }

  cord<-GetTissueCoordinates(obj)
  nn.idx<-RANN::nn2(cord[,1:2],k=k_spatial)$nn.idx
  if (usepca){exp<-t(Embeddings(obj, reduction = "pca"))}else{
    if (is_v5) {
      exp <- GetAssayData(obj, assay = assay, layer = "data")
    } else {
      exp <- GetAssayData(obj, assay = assay, slot = "data")
    }}

  # Number of cells
  num_cells <- ncol(exp)
  num_neighbors <- k_spatial - 1  # Exclude self

  # Create a sparse matrix for storing the neighbor information
  # Each row of nn.idx contains the indices of neighbors for the corresponding cell
  # So, we will create a sparse matrix where each entry corresponds to a neighbor of a cell
  # This matrix will have 1s for valid neighbor positions (we'll use it as a "weight" matrix)
  nn_matrix <- Matrix::sparseMatrix(j = rep(1:num_cells, times=num_neighbors),
                                    i = as.vector(nn.idx[, -1]),
                                    x = 1/num_neighbors, dims = c(num_cells, num_cells))

  # Multiply the expression matrix (exp) with the nn_matrix to get the sum of neighbor expression values
  nn.exp<- exp %*% (nn_matrix)
  message("Complete step 1: Generate the microenvironment profile")

  # Use RANN's nn2 directly and avoid repeated transpositions

  pcaembed <- Embeddings(obj, reduction = "pca")


  nnmtx <- RANN::nn2(pcaembed, k = 2 * k_nn)







  res<-find_knn_rcpp(as.matrix(nn.exp),nnmtx$nn.idx)
  index_reorder <- res$indices[, 1:k_nn]
  message("Complete step 2: Refine the expression-NN based on the microenvironment profile")


  rownum<-nrow(index_reorder)
  ##adaptive
  if (K_adp==T) {


    ###optimized version####

    nndist_reorder<-res$distances[,1:(k_nn+1)]

    denom   <- k_nn - 1 - delta
    thresh  <- (rowSums( sqrt(nndist_reorder) ) / denom)^2   # length = nrow
    ind     <- sweep(nndist_reorder, 1, thresh, FUN = "<")
    ind<-ind[,-ncol(ind)]

    rowind<-rep(index_reorder[,1],apply(ind,1,sum))
    colind<-t(index_reorder)[t(ind)]
  }else{

    rowind<-rep(index_reorder[,1],each=k_nn)
    colind<-t(index_reorder)

  }

  tmpsparse<- Matrix::sparseMatrix(
    i = rowind,
    j = colind,
    x = rep(1,length(rowind)),
    dims = c(rownum,rownum),
    dimnames = list(colnames(obj),colnames(obj)))

  snnmatrix<-tmpsparse %*% Matrix::t(tmpsparse)
  snnmatrix@x<-snnmatrix@x/(k_nn*2-snnmatrix@x)
  snnmatrix<-as(snnmatrix,"TsparseMatrix")


  keep<-snnmatrix@x>prune.SNN
  snnmatrix@i<-snnmatrix@i[keep]
  snnmatrix@j<-snnmatrix@j[keep]
  snnmatrix@x<-snnmatrix@x[keep]
  snn.graph<-as.Graph(snnmatrix)

  DefaultAssay(snn.graph) <- assay
  obj@graphs$MEcell<-snn.graph
  return(obj)

}


