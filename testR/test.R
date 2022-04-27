rm(list = ls())
load_all()
library(tidyverse)
library(Seurat)

reference <- readRDS("~/github/seurat/data/NKT_ref_batches_2022-04-20.rds")
# reference <- NormalizeData(reference)
#
query <- readRDS("~/github/seurat/data/NKT_query_batch_2022-04-20.rds")
# query <- NormalizeData(query)
#
#
DEGsT <- read_tsv("~/github/seurat/data/snn-single-markers.tsv")
DEGs <- unique(DEGsT %>% pull(gene))
#
# A.L <- FindTransferAnchors_MB(reference = reference, query = query, features = DEGs, reduction = "pcaproject")
# saveRDS(A.L, "~/github/seurat/data/A.L")

A.L <- readRDS("~/github/seurat/data/A.L")


# gamma.L <- GetRefBatchWeight_MB(A.L, reference@meta.data$seurat_clusters)
# saveRDS(gamma.L, "~/github/seurat/data/gamma.L")
gamma.L <- readRDS("~/github/seurat/data/gamma.L")

query <- TransferData_MB(
  anchorsetL = A.L[gamma.L$batch[1:20]],
  refdata = "seurat_clusters",
  k.weight = 30,
  reference = reference,
  query = query )

# sort(unlist(gamma.L), decreasing = T)
