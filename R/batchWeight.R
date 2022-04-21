#' Find transfer anchors list
#'
#' Find a set of anchors between each reference batch and query object. These
#' anchors can later be used to transfer data from the reference to
#' query object using the \code{\link{TransferData}} object.
#'
#' @param reference
#' @param query
#' @param batch
#' @param normalization.method
#' @param recompute.residuals
#' @param reference.assay
#' @param reference.neighbors
#' @param query.assay
#' @param reduction
#' @param reference.reduction
#' @param project.query
#' @param features
#' @param scale
#' @param npcs
#' @param l2.norm
#' @param dims
#' @param k.anchor
#' @param k.filter
#' @param k.score
#' @param max.features
#' @param nn.method
#' @param n.trees
#' @param eps
#' @param approx.pca
#' @param mapping.score.k
#' @param verbose
#'
#' @return
#' @export
#'
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
