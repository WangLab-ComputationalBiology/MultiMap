ValidateParams_TransferData <- function(
  anchorset,
  combined.ob,
  anchors,
  reference.cells,
  query.cells,
  reference,
  query,
  refdata,
  weight.reduction,
  l2.norm,
  dims,
  k.weight,
  sd.weight,
  eps,
  n.trees,
  verbose,
  slot,
  prediction.assay,
  label.transfer) {
  if (!inherits(x = refdata, what = "list")) {
    refdata <- list(id = refdata)
  }
  for (i in 1:length(x = refdata)) {
    if (inherits(x = refdata[[i]], what = c("character", "factor"))) {
      # check is it's in the reference object
      if (length(x = refdata[[i]]) == 1) {
        if (is.null(x = reference)) {
          warning("If providing a single string to refdata element number ", i,
                  ", please provide the reference object. Skipping element ", i,
                  ".", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
          next
        }
        if (refdata[[i]] %in% Assays(object = reference)) {
          refdata[[i]] <- GetAssayData(object = reference, assay = refdata[[i]])
          colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
          label.transfer[[i]] <- FALSE
          next
        } else if (refdata[[i]] %in% colnames(x = reference[[]])) {
          refdata[[i]] <- reference[[refdata[[i]]]][, 1]
        } else {
          warning("Element number ", i, " provided to refdata does not exist in ",
                  "the provided reference object.", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
          next
        }
      } else if (length(x = refdata[[i]]) != length(x = reference.cells)) {
        warning("Please provide a vector that is the same length as the number ",
                "of reference cells used in anchor finding.\n",
                "Length of vector provided: ", length(x = refdata[[i]]), "\n",
                "Length of vector required: ", length(x = reference.cells),
                "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
        refdata[[i]] <- FALSE
      }
      label.transfer[[i]] <- TRUE
    } else if (inherits(x = refdata[[i]], what = c("dgCMatrix", "matrix"))) {
      if (ncol(x = refdata[[i]]) != length(x = reference.cells)) {
        warning("Please provide a matrix that has the same number of columns as ",
                "the number of reference cells used in anchor finding.\n",
                "Number of columns in provided matrix : ", ncol(x = refdata[[i]]), "\n",
                "Number of columns required           : ", length(x = reference.cells),
                "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
        refdata[[i]] <- FALSE
      } else {
        colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
        if (any(!colnames(x = refdata[[i]]) == reference.cells)) {
          if (any(!colnames(x = refdata[[i]]) %in% reference.cells) || any(!reference.cells %in% colnames(x = refdata[[i]]))) {
            warning("Some (or all) of the column names of the provided refdata ",
                    "don't match the reference cells used in anchor finding ",
                    "\nSkipping element", i, ".", call. = FALSE, immediate. = TRUE)
            refdata[[i]] <- FALSE
          } else {
            refdata[[i]] <- refdata[[i]][, reference.cells]
          }
        }
      }
      if (!slot %in% c("counts", "data")) {
        stop("Please specify slot as either 'counts' or 'data'.")
      }
      label.transfer[[i]] <- FALSE
    } else {
      warning("Please provide either a vector (character or factor) for label ",
              "transfer or a matrix for feature transfer. \nType provided: ",
              class(x = refdata[[i]]))
      refdata[[i]] <- FALSE
    }
    if (names(x = refdata)[i] == "") {
      possible.names <- make.unique(names = c(names(x = refdata), paste0("e", i)))
      names(x = refdata)[i] <- possible.names[length(x = possible.names)]
      if (verbose) {
        message("refdata element ", i, " is not named. Setting name as ", names(x = refdata)[i])
      }
    }
  }
  ModifyParam(param = "label.transfer", value = label.transfer)
  if (all(unlist(x = lapply(X = refdata, FUN = isFALSE)))) {
    stop("None of the provided refdata elements are valid.", call. = FALSE)
  }
  ModifyParam(param = "refdata", value = refdata)
  valid.weight.reduction <- c("pcaproject", "pca", "cca", "rpca.ref","lsiproject", "lsi")
  if (!inherits(x = weight.reduction, "DimReduc")) {
    if (!weight.reduction %in% valid.weight.reduction) {
      stop("Please provide one of ", paste(valid.weight.reduction, collapse = ", "), " or a custom DimReduc to ",
           "the weight.reduction parameter.", call. = FALSE)
    }
    if (weight.reduction %in% c("pcaproject", "cca", "rpca.ref", "lsiproject") &&
        !weight.reduction %in% Reductions(object = combined.ob)) {
      stop("Specified weight.reduction (", weight.reduction, ") is not present ",
           "in the provided anchorset.", call. = FALSE)
    }
    if (weight.reduction %in% c("pca", "lsi") && is.null(x = query)) {
      stop("To use an internal PCA on the query only for weight.reduction, ",
           "please provide the query object.", call. = FALSE)
    }
  }
  if (inherits(x = weight.reduction, "DimReduc")) {
    if (is.null(x = dims)) {
      stop("Please specify dims", call. = FALSE)
    }
    if (max(dims) > ncol(x = weight.reduction)) {
      stop("The max of dims specified (", max(dims), ") is greater than the ",
           "number of dimensions in the given DimReduc (",
           ncol(x = weight.reduction), ").", call. = FALSE)
    }
  } else {
    if (is.null(x = dims)) {
      ModifyParam(param = "dims", value = 1:length(x = slot(object = anchorset, name = "command")$dims))
    }
  }

  if (!is.null(x = query)) {
    if (!isTRUE(x = all.equal(
      target = gsub(pattern = "_query", replacement = "", x = query.cells),
      current = Cells(x = query),
      check.attributes = FALSE)
      )) {
      stop("Query object provided contains a different set of cells from the ",
           "query used to construct the AnchorSet provided.", call. = FALSE)
    }
  }
  if(k.weight > nrow(x = anchors)) {
    stop("Please set k.weight to be smaller than the number of anchors (",
         nrow(x = anchors), ").", call. = FALSE)
  }
}

# here apply pca on query data, move to
SetObj_SB <- function(anchorset,
                      refdata,
                      reference.total = NULL,
                      query = NULL,
                      weight.reduction = 'pcaproject',
                      l2.norm = FALSE,
                      dims = NULL,
                      k.weight = 50,
                      sd.weight = 1,
                      eps = 0,
                      n.trees = 50,
                      verbose = TRUE,
                      slot = "data",
                      prediction.assay = FALSE,
                      store.weights = TRUE
                      ) {
    # step1 : filter out good batches #########################################
    # Get batch score/weight  #################################################
    # outside test ############################################################
    # step2 : generate weights ################################################
    # pre-process #############################################################
    # pp step 1: parameter validation #########################################
    # pp step 2: data format correction #######################################
    # this step need to nn2 steps #############################################
    # nn2 step1: for each batch find neighbors for normalization ##############
    # nn2 step2: for all anchors weight normalization #########################
    # step3 : transfer labels #################################################
    # this step need no integration operations ################################


    Idents(reference.total) <- reference.total$batch
    reference <- subset(reference.total, idents = names(anchorset))

    combined.ob <- slot(object = anchorset, name = "object.list")[[1]]
    anchors <- slot(object = anchorset, name = "anchors")

    # all reference cells #####################################################

    reference.cells <- slot(object = anchorset, name = "reference.cells")

    # all query cells #########################################################

    query.cells <- slot(object = anchorset, name = "query.cells")
    label.transfer <- list()

    # here needs a new ValidateParams func for batches transfer ###############

    ValidateParams_TransferData(
        anchorset = anchorset,
        combined.ob = combined.ob,
        anchors = anchors,
        reference.cells = reference.cells,
        query.cells = query.cells,
        refdata = refdata,
        reference = reference,
        query = query,
        weight.reduction = weight.reduction,
        l2.norm = l2.norm,
        dims = dims,
        k.weight = k.weight,
        sd.weight = sd.weight,
        eps = eps,
        n.trees = n.trees,
        verbose = verbose,
        slot = slot,
        prediction.assay = prediction.assay,
        label.transfer = label.transfer)


    combined.ob <- SetIntegrationData(
        object = combined.ob,
        integration.name = "integrated",
        slot = 'anchors',
        new.data = anchors)

    combined.ob <- SetIntegrationData(
        object = combined.ob,
        integration.name = "integrated",
        slot = 'neighbors',
        new.data = list('cells1' = reference.cells, 'cells2' = query.cells))

    return(combined.ob)
}



SetObj_MB <- function(anchorsetL,
                      refdata,
                      reference = NULL,
                      query = NULL,
                      weight.reduction = 'pcaproject',
                      l2.norm = FALSE,
                      dims = NULL,
                      k.weight = 50,
                      sd.weight = 1,
                      eps = 0,
                      n.trees = 50,
                      verbose = TRUE,
                      slot = "data",
                      prediction.assay = FALSE,
                      store.weights = TRUE
                      ){

    combined.ob.L <- lapply(anchorsetL,
                            SetObj_SB,
                            refdata,
                            reference.total = reference,
                            query,
                            weight.reduction,
                            l2.norm,
                            dims,
                            k.weight,
                            sd.weight,
                            eps,
                            n.trees,
                            verbose,
                            slot,
                            prediction.assay,
                            store.weights)
   return(combined.ob.L)
}


FindDistances_SB <- function(object,
                             reduction = NULL,
                             assay = NULL,
                             integration.name = 'integrated',
                             dims = 1:10,
                             features = NULL,
                             k = 300,
                             sd.weight = 1,
                             nn.method = "annoy",
                             n.trees = 50,
                             eps = 0,
                             verbose = TRUE){
  if (verbose) {
    message("Finding integration vector weights")
  }
  if (is.null(x = reduction) & is.null(x = features)) {
    stop("Need to specify either dimension reduction object or a set of features")
  }
  assay <- assay %||% DefaultAssay(object = object)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2

    message("get neighbors")
  anchors <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors'
  )
    message("get anchors")
    anchors.cells2 <- unique(x = nn.cells2[anchors[, "cell2"]])
    if (is.null(x = features)) {
      data.use <- Embeddings(reduction)[nn.cells2, dims]
    } else {
      data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, nn.cells2])
    }

    message("begin NNhelper")
    knn_2_2 <- NNHelper(
      data = data.use[anchors.cells2, ],
      query = data.use,
      k = k,
      method = nn.method,
      n.trees = n.trees,
      eps = eps
    )
    message("get NNhelper")

  distances <- Distances(object = knn_2_2)
  distances <- 1 - (distances / distances[, ncol(x = distances)])
  cell.index <- Indices(object = knn_2_2)

  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'distances',
    new.data = distances
  )

  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'cell.index',
    new.data = cell.index
  )

  return(object)
}



FindDistances_MB <- function(object.L,
                             reduction = NULL,
                             assay = NULL,
                             integration.name = 'integrated',
                             dims = 1:10,
                             features = NULL,
                             k = 300,
                             sd.weight = 1,
                             nn.method = "annoy",
                             n.trees = 50,
                             eps = 0,
                             verbose = TRUE
                             ){
    combined.ob.L <- lapply(object.L,
                           FindDistances_SB,
                           reduction,
                           assay,
                           integration.name,
                           dims,
                           features,
                           k,
                           sd.weight,
                           nn.method,
                           n.trees,
                           eps,
                           verbose)
    return(combined.ob.L)
}

GetComAnchorsC2_MB <- function(combined.ob.L, integration.name){
    C2L <- lapply(combined.ob.L, GetComAnchorsC2_SB, integration.name)
    return(unlist(C2L))
}

GetComAnchorsC2_SB <- function(object, integration.name){
    anchors <- GetIntegrationData(
        object = object,
        integration.name = integration.name,
        slot = 'anchors'
    )
    return(anchors[, "cell2"])
}

GetComAnchors_MB <- function(combined.ob.L, integration.name){
    CL <- lapply(combined.ob.L, GetComAnchors_SB, integration.name)
    anchors <- bind_rows(CL, .id = "column_label")
    return(anchors)
}

GetComAnchors_SB <- function(object, integration.name){
    anchors <- GetIntegrationData(
        object = object,
        integration.name = integration.name,
        slot = 'anchors'
    )
    return(as.data.frame(anchors))
}

getIntegrationRownames_MB <- function(combined.ob.L, integration.name){
    rownamesL <- lapply(combined.ob.L, getIntegrationRownames_SB, integration.name)
    return(unlist(rownamesL))
}

getIntegrationRownames_SB <- function(object, integration.name){
    integration.matrix <- GetIntegrationData(
        object = object,
        integration.name = integration.name,
        slot = "integration.matrix")
    return(rownames(integration.matrix))
}

getAnchorScore_MB <- function(combined.ob.L, integration.name){
    S2L <- lapply(combined.ob.L, getAnchorScore_SB, integration.name)
    return(unlist(S2L))
}
getAnchorScore_SB <- function(object, integration.name){
    anchors <- GetIntegrationData(
        object = object,
        integration.name = integration.name,
        slot = 'anchors'
    )
    return(anchors[, "score"])
}

FindIntegrationMatrix <- function(
  object,
  assay = NULL,
  integration.name = 'integrated',
  features.integrate = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  anchors <- GetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors'
  )
  if (verbose) {
    message("Finding integration vectors")
  }

  features.integrate <- features.integrate %||% rownames(
    x = GetAssayData(object = object, assay = assay, slot = "data")
  )
  data.use1 <- t(x = as.matrix(GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells1]
  ))
  data.use2 <- t(x = as.matrix(GetAssayData(
    object = object,
    assay = assay,
    slot = "data")[features.integrate, nn.cells2]
  ))
  anchors1 <- nn.cells1[anchors[, "cell1"]]
  anchors2 <- nn.cells2[anchors[, "cell2"]]
  data.use1 <- data.use1[anchors1, ]
  data.use2 <- data.use2[anchors2, ]
  integration.matrix <- data.use2 - data.use1
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'integration.matrix',
    new.data = integration.matrix
  )
  return(object)
}

FindWeights_MB <- function(combined.ob.L,
                           query,
                           refdata,
                           reduction = NULL,
                           assay = NULL,
                           integration.name = 'integrated',
                           dims = 1:10,
                           features = NULL,
                           k = 300,
                           sd.weight = 1,
                           nn.method = "annoy",
                           n.trees = 50,
                           eps = 0,
                           reverse = FALSE,
                           verbose = TRUE) {
    # input: query cells, all anchors' distance list, all anchors' scoer list, etc ############################
    # step1: get distance  ####################################################
    object = combined.ob.L[[1]]
    if (verbose) {
        message("Finding integration vector weights")
    }
    message("find features")
    if (is.null(x = reduction) & is.null(x = features)) {
        stop("Need to specify either dimension reduction object or a set of features")
    }
    message("end find features")
    assay <- assay %||% DefaultAssay(object = object)
    neighbors <- GetIntegrationData(object = object,
                                    integration.name = integration.name,
                                    slot = 'neighbors')
    nn.cells1 <- neighbors$cells1
    nn.cells2 <- neighbors$cells2

    anchors.cells2 <- unique(x = nn.cells2[GetComAnchorsC2_MB(combined.ob.L, integration.name)])

    if (is.null(x = features)) {
        data.use <- Embeddings(reduction)[nn.cells2, dims]
    } else {
        data.use <- t(x = GetAssayData(object = object,
                                       slot = 'data',
                                       assay = assay)[features, nn.cells2])
    }

    knn_2_2 <- NNHelper(
        data = data.use[anchors.cells2, ],
        query = data.use,
        k = k,
        method = nn.method,
        n.trees = n.trees,
        eps = eps)

    distances <- Distances(object = knn_2_2)
    distances <- 1 - (distances / distances[, ncol(x = distances)])
    cell.index <- Indices(object = knn_2_2)

    # step 2 get integration matrix  ##########################################
    combined.ob.L <- lapply(combined.ob.L, FindIntegrationMatrix,
                            assay = "RNA", integration.name = 'integrated',
                            features.integrate = NULL, verbose = TRUE)

    integration_matrix_rownames_MB <- getIntegrationRownames_MB(combined.ob.L, integration.name = integration.name)
    anchor_score_MB <- getAnchorScore_MB(combined.ob.L, integration.name = integration.name)

    weights <- FindWeightsC(
        # query cells #########################################################
        cells2 = 0:(length(x = nn.cells2) - 1),
        # distance list, for each batch #######################################
        distances = as.matrix(x = distances),
        # distance list, for each batch #######################################
        anchor_cells2 = anchors.cells2,
        # return integration_matrix_rownames.size * cell2.size matrix #########
        # for label prediction ################################################
        integration_matrix_rownames = integration_matrix_rownames_MB,
        cell_index = cell.index,
        # anchor_score input as list ##########################################
        # the same order and size of integration_matrix_rownames ##############
        anchor_score = anchor_score_MB,
        min_dist = 0,
        sd = sd.weight,
        display_progress = verbose)

    rownames(weights) <- integration_matrix_rownames_MB
    colnames(weights) <- nn.cells2

    # prediction: step1  get nn.cell1 which <=> nn.cells2
    nn.anchors <- GetComAnchors_MB(combined.ob.L, integration.name)

    query.cells <- unname(obj = sapply(
      X = nn.cells2,
      FUN = function(x) gsub(pattern = "_query", replacement = "", x = x)
    ))

    # step2: get annotation of nn.cell1 from refdata
    nn.anchors$id1 <- refdata[[1]][nn.anchors[,"cell1"]]

    # step3: generate prediction matrix
    reference.ids <- factor(x = nn.anchors$id1, levels = unique(x = refdata[[1]]))
    possible.ids <- levels(x = reference.ids)

    prediction.mat <- matrix(nrow = nrow(x = nn.anchors), ncol = length(x = possible.ids), data = 0)
    for(i in 1:length(x = possible.ids)){
      prediction.mat[which(reference.ids == possible.ids[i]), i] = 1
    }

    if(verbose){
      message("Predicting cell labels")
    }


    # step4: get prediction.scores

    browser()
    prediction.scores <- Matrix::t(x = weights) %*% prediction.mat

    # step5: generate metadata

    colnames(x = prediction.scores) <- possible.ids
    rownames(x = prediction.scores) <- query.cells
    prediction.ids <- possible.ids[apply(X = prediction.scores, MARGIN = 1, FUN = which.max)]
    prediction.ids <- as.character(prediction.ids)
    prediction.max <- apply(X = prediction.scores, MARGIN = 1, FUN = max)
    if (is.null(x = query)){
      prediction.scores <- cbind(prediction.scores, max = prediction.max)
    }
    predictions <- data.frame(
      predicted.id = prediction.ids,
      prediction.score = as.matrix(prediction.scores),
      row.names = query.cells,
      stringsAsFactors = FALSE
    )

    query <- AddMetaData(object = query, metadata = prediction.max, col.name = paste0("predicted.score"))
    query <- AddMetaData(object = query, metadata = prediction.ids, col.name = paste0("predicted"))

    return(query)
}



TransferData_MB <- function(# input anchor set list  #################
                            anchorsetL,
                            refdata,
                            reference = NULL,
                            query = NULL,
                            features = NULL,
                            assay = "RNA",
                            integration.name = 'integrated',
                            weight.reduction = 'pcaproject',
                            l2.norm = FALSE,
                            dims = 1:10,
                            nn.method = "annoy",
                            k.weight = 50,
                            sd.weight = 1,
                            eps = 0,
                            n.trees = 50,
                            verbose = TRUE,
                            slot = "data",
                            prediction.assay = FALSE,
                            store.weights = TRUE
                            ) {


    # step1  parameter validation ###################################################################
    message("step1")
    combined.ob.L <- SetObj_MB(anchorsetL,
                               refdata,
                               reference,
                               query,
                               weight.reduction,
                               l2.norm,
                               dims,
                               k.weight,
                               sd.weight,
                               eps,
                               n.trees,
                               verbose,
                               slot,
                               prediction.assay,
                               store.weights)


  if (!inherits(x = refdata, what = "list")) {
    refdata <- list(id = refdata)
  }
  for (i in 1:length(x = refdata)) {
    if (inherits(x = refdata[[i]], what = c("character", "factor"))) {
      # check is it's in the reference object
      if (length(x = refdata[[i]]) == 1) {
        if (is.null(x = reference)) {
          warning("If providing a single string to refdata element number ", i,
                  ", please provide the reference object. Skipping element ", i,
                  ".", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
          next
        }
        if (refdata[[i]] %in% Assays(object = reference)) {
          refdata[[i]] <- GetAssayData(object = reference, assay = refdata[[i]])
          colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
          next
        } else if (refdata[[i]] %in% colnames(x = reference[[]])) {
          refdata[[i]] <- reference[[refdata[[i]]]][, 1]
        } else {
          warning("Element number ", i, " provided to refdata does not exist in ",
                  "the provided reference object.", call. = FALSE, immediate. = TRUE)
          refdata[[i]] <- FALSE
          next
        }
      } else if (length(x = refdata[[i]]) != length(x = reference.cells)) {
        warning("Please provide a vector that is the same length as the number ",
                "of reference cells used in anchor finding.\n",
                "Length of vector provided: ", length(x = refdata[[i]]), "\n",
                "Length of vector required: ", length(x = reference.cells),
                "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
        refdata[[i]] <- FALSE
      }
    } else if (inherits(x = refdata[[i]], what = c("dgCMatrix", "matrix"))) {
      if (ncol(x = refdata[[i]]) != length(x = reference.cells)) {
        warning("Please provide a matrix that has the same number of columns as ",
                "the number of reference cells used in anchor finding.\n",
                "Number of columns in provided matrix : ", ncol(x = refdata[[i]]), "\n",
                "Number of columns required           : ", length(x = reference.cells),
                "\nSkipping element ", i, ".", call. = FALSE, immediate. = TRUE)
        refdata[[i]] <- FALSE
      } else {
        colnames(x = refdata[[i]]) <- paste0(colnames(x = refdata[[i]]), "_reference")
        if (any(!colnames(x = refdata[[i]]) == reference.cells)) {
          if (any(!colnames(x = refdata[[i]]) %in% reference.cells) || any(!reference.cells %in% colnames(x = refdata[[i]]))) {
            warning("Some (or all) of the column names of the provided refdata ",
                    "don't match the reference cells used in anchor finding ",
                    "\nSkipping element", i, ".", call. = FALSE, immediate. = TRUE)
            refdata[[i]] <- FALSE
          } else {
            refdata[[i]] <- refdata[[i]][, reference.cells]
          }
        }
      }
      if (!slot %in% c("counts", "data")) {
        stop("Please specify slot as either 'counts' or 'data'.")
      }
    } else {
      warning("Please provide either a vector (character or factor) for label ",
              "transfer or a matrix for feature transfer. \nType provided: ",
              class(x = refdata[[i]]))
      refdata[[i]] <- FALSE
    }
    if (names(x = refdata)[i] == "") {
      possible.names <- make.unique(names = c(names(x = refdata), paste0("e", i)))
      names(x = refdata)[i] <- possible.names[length(x = possible.names)]
      if (verbose) {
        message("refdata element ", i, " is not named. Setting name as ", names(x = refdata)[i])
      }
    }
  }
    # step2 ###################################################################
    # step2.1 apply pca on query data, then pass it as parameter to
    if (!inherits(x = weight.reduction, what = "DimReduc") && weight.reduction == 'pca') {
        if (verbose) {
            message("Running PCA on query dataset")
        }
        features <- slot(object = anchorset, name = "anchor.features")
        query.ob <- query
        query.ob <- ScaleData(object = query.ob, features = features, verbose = FALSE)
        query.ob <- RunPCA(object = query.ob, npcs = max(dims), features = features, verbose = FALSE)
        rm(query.ob)
    }

    if (!inherits(x = weight.reduction, what = "DimReduc") && weight.reduction == "lsi") {
        if (!("lsi" %in% Reductions(object = query))) {
            stop("Requested lsi for weight.reduction, but lsi not stored in query object.")
        } else {
            weight.reduction <- query[["lsi"]]
        }
    }

    if (inherits(x = weight.reduction, what = "DimReduc")) {
        weight.reduction <- RenameCells(
            object = weight.reduction,
            new.names = paste0(Cells(x = weight.reduction), "_query"))
    } else {
        if (l2.norm) {
            weight.reduction.l2 <- paste0(weight.reduction, ".l2")
            if (weight.reduction.l2 %in% Reductions(object = combined.ob.L[[1]])) {
                combined.ob.L[[1]] <- L2Dim(object = combined.ob.L[[1]], reduction = weight.reduction)
            }
            weight.reduction <- weight.reduction.l2
        }
        weight.reduction <- combined.ob.L[[1]][[weight.reduction]]
    }

    if (max(dims) > ncol(x = weight.reduction)) {
        stop("dims is larger than the number of available dimensions in ",
             "weight.reduction (", ncol(x = weight.reduction), ").", call. = FALSE)
    }

    # step 2: get distance
    # message("step2")
    # combined.ob.L <- FindDistances_MB(combined.ob.L,
    #                                   reduction = weight.reduction,
    #                                   assay = assay,
    #                                   integration.name = 'integrated',
    #                                   dims,
    #                                   features = features,
    #                                   k = k.weight,
    #                                   sd.weight,
    #                                   nn.method,
    #                                   n.trees,
    #                                   eps,
    #                                   verbose)
#
    # # step3 ###################################################################
#
    message("step3")
    query <- FindWeights_MB(combined.ob.L,
                            query = query,
                              refdata = refdata,
                                    reduction = weight.reduction,
                                    assay = assay,
                                    integration.name = integration.name,
                                    dims = dims,
                                    features = features,
                                    k = k.weight,
                                    sd.weight = sd.weight,
                                    nn.method = nn.method,
                                    n.trees = n.trees,
                                    eps = eps,
                                    verbose = verbose)
#
#
    return(query)
}
