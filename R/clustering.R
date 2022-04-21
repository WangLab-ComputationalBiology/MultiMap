# Internal helper function to dispatch to various neighbor finding methods
#
# @param data Input data
# @param query Data to query against data
# @param k Number of nearest neighbors to compute
# @param method Nearest neighbor method to use: "rann", "annoy"
# @param cache.index Store algorithm index with results for reuse
# @param ... additional parameters to specific neighbor finding method
#
#' @importFrom methods new
#' @importClassesFrom SeuratObject Neighbor
#
NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
        do.call(what = 'AnnoyNN', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'annoy'")
    )
  )
  n.ob <- new(
    Class = 'Neighbor',
    nn.idx = results$nn.idx,
    nn.dist = results$nn.dists,
    alg.info = results$alg.info %||% list(),
    cell.names = rownames(x = query)
  )
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }
  return(n.ob)
}
