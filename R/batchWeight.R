#' Find transfer anchors list
#'
#' Find a set of anchors between each reference batch and query object. These
#' anchors can later be used to transfer data from the reference to
#' query object using the \code{\link{TransferData}} object.
#'
#' @param reference \code{\link{Seurat}} object to use as the reference
#' @param query \code{\link{Seurat}} object to use as the query
#' @param batch batch column in meta.data of Seurat object
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT.
#' @param recompute.residuals If using SCT as a normalization method, compute
#' query Pearson residuals using the reference SCT model parameters.
#' @param reference.assay Name of the Assay to use from reference
#' @param reference.neighbors Name of the Neighbor to use from the reference.
#' Optionally enables reuse of precomputed neighbors.
#' @param query.assay Name of the Assay to use from query
#' @param reduction Dimensional reduction to perform when finding anchors.
#' Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the reference onto the query. We
#'    recommend using PCA when reference and query datasets are from scRNA-seq}
#'    \item{lsiproject: Project the LSI from the reference onto the query. We
#'    recommend using LSI when reference and query datasets are from scATAC-seq.
#'    This requires that LSI has been computed for the reference dataset, and the
#'    same features (eg, peaks or genome bins) are present in both the reference
#'    and query. See \code{\link[Signac]{RunTFIDF}} and
#'    \code{\link[Signac]{RunSVD}}}
#'    \item{rpca: Project the PCA from the reference onto the query, and the PCA
#'    from the query onto the reference (reciprocal PCA projection).}
#'    \item{cca: Run a CCA on the reference and query }
#' }
#' @param reference.reduction Name of dimensional reduction to use from the
#' reference if running the pcaproject workflow. Optionally enables reuse of
#' precomputed reference dimensional reduction. If NULL (default), use a PCA
#' computed on the reference object.
#' @param project.query Project the PCA from the query dataset onto the
#' reference. Use only in rare cases where the query dataset has a much larger
#' cell number, but the reference dataset has a unique assay for transfer. In
#' this case, the default features will be set to the variable features of the
#' query object that are alos present in the reference.
#' @param features Features to use for dimensional reduction. If not specified,
#' set as variable features of the reference object which are also present in
#' the query.
#' @param scale Scale query data.
#' @param npcs Number of PCs to compute on reference if reference.reduction is
#' not provided.
#' @param l2.norm Perform L2 normalization on the cell embeddings after
#' dimensional reduction
#' @param dims Which dimensions to use from the reduction to specify the
#' neighbor search space
#' @param k.anchor How many neighbors (k) to use when finding anchors
#' @param k.filter How many neighbors (k) to use when filtering anchors. Set to
#' NA to turn off filtering.
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param max.features The maximum number of features to use when specifying the
#' neighborhood search space in the anchor filtering
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param eps Error bound on the neighbor finding algorithm (from
#' \code{\link{RANN}} or \code{\link{RcppAnnoy}})
#' @param approx.pca Use truncated singular value decomposition to approximate
#' PCA
#' @param mapping.score.k Compute and store nearest k query neighbors in the
#' AnchorSet object that is returned. You can optionally set this if you plan
#' on computing the mapping score and want to enable reuse of some downstream
#' neighbor calculations to make the mapping score function more efficient.
#' @param verbose Print progress bars and output
#'
#' @return Returns an \code{AnchorSet} object list.
#'
#' @export
#' @examples
FindTransferAnchors_MB <- function(reference,
                                   query,
                                   batch = "batch",
                                   normalization.method = "LogNormalize",
                                   recompute.residuals = TRUE,
                                   reference.assay = NULL,
                                   reference.neighbors = NULL,
                                   query.assay = NULL,
                                   reduction = "pcaproject",
                                   reference.reduction = NULL,
                                   project.query = FALSE,
                                   features = NULL,
                                   scale = TRUE,
                                   npcs = 30,
                                   l2.norm = TRUE,
                                   dims = 1:30,
                                   k.anchor = 5,
                                   k.filter = 200,
                                   k.score = 30,
                                   max.features = 200,
                                   nn.method = "annoy",
                                   n.trees = 50,
                                   eps = 0,
                                   approx.pca = TRUE,
                                   mapping.score.k = NULL,
                                   verbose = TRUE ) {
  anchorsL <- list()
  sortedBatch <- sort(table(reference@meta.data[,batch]), decreasing = T)
  for(i in 1:length(sortedBatch)){
    RefName <- names(sortedBatch)[i]
    tempRefObj <- subset(reference,
                         cells = Cells(reference)[reference@meta.data[,batch] == RefName])
    tryCatch({
      temp.anchors <- FindTransferAnchors(reference = tempRefObj,
                                          query = query,
                                          normalization.method = normalization.method,
                                          recompute.residuals = recompute.residuals,
                                          reference.assay = reference.assay,
                                          reference.neighbors = reference.neighbors,
                                          query.assay = query.assay,
                                          reduction = reduction,
                                          reference.reduction = reference.reduction,
                                          project.query = project.query,
                                          features = features,
                                          scale = scale,
                                          npcs = npcs,
                                          l2.norm = l2.norm,
                                          dims = dims,
                                          k.anchor = k.anchor,
                                          k.filter = k.filter,
                                          k.score = k.score,
                                          max.features = max.features,
                                          nn.method = nn.method,
                                          n.trees = n.trees,
                                          eps = eps,
                                          approx.pca = approx.pca,
                                          mapping.score.k = mapping.score.k,
                                          verbose = verbose)
      anchorsL[[RefName]] <- temp.anchors
    },
    error = function(e){print(e)}
    )
  }
  return(anchorsL)
}



#' Get reference quality score gamma for each batch
#'
#' @param anchorsetL an \code{AnchorSet} object list.
#' @param celltype cell type vector for each cell
#' @param beta beta parameter
#' @param alpha alpha parameter
#'
#' @return gamma score list
#' @export
#'
#' @examples
GetRefBatchWeight_MB <- function(anchorsetL,
                              celltype,
                              beta = 0.5,
                              alpha = 0.5){
    ## cell1 cell2     score
    ## [1,]     7   766 0.6923077
    gammaL <- lapply(anchorsetL, GetBatchWeight, celltype, beta)
    gammaT <- bind_rows(gammaL) %>%
      mutate(batch = names(gammaL)) %>%
      mutate(Snorm = S/max(S)) %>%
      mutate(ICnorm = IC/max(IC)) %>%
      mutate(Gamma = alpha * log(ICnorm) + (1 - alpha) * log(Snorm)) %>%
      arrange(desc(Gamma)) %>%
      select(batch, IC, ICnorm, S, Snorm, Gamma)

    return(gammaT)
}

#' Get gamma for a batch
#'
#' @param anchorset an \code{AnchorSet} object
#' @param celltype cell type vector for each cell
#' @param beta beta parameter
#' @param alpha alpha parameter
#'
#' @return gamma score
#'
#' @examples
GetBatchWeight <- function(anchorset, celltype, beta){
    scoreT <- as.tibble(anchorset@anchors)
    K <- length(unique(x = celltype))
    scoreT$celltype <- celltype[scoreT$cell1]

    scoreT <- scoreT %>%
        mutate(goodscore = (score > beta) * score)
    TotalGoodScore <- sum(scoreT$goodscore)

    IC <- scoreT %>%
      group_by(celltype) %>%
      summarise(GS = sum(goodscore)) %>%
      ungroup %>%
      mutate(CIC = - GS / TotalGoodScore * log(GS / TotalGoodScore)) %>%
      replace_na(list(CIC = 0)) %>%
      summarise(IC = sum(CIC) / K) %>%
      pull(IC)

    S <- TotalGoodScore / K

    return(tibble(IC=IC, S=S))
}

