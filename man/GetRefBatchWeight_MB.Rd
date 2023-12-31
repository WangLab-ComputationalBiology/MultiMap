% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batchWeight.R
\name{GetRefBatchWeight_MB}
\alias{GetRefBatchWeight_MB}
\title{Get reference quality score gamma for each batch}
\usage{
GetRefBatchWeight_MB(anchorsetL, celltype, beta = 0.5, alpha = 0.5)
}
\arguments{
\item{anchorsetL}{an \code{AnchorSet} object list.}

\item{celltype}{cell type vector for each cell}

\item{beta}{beta parameter}

\item{alpha}{alpha parameter}
}
\value{
gamma score table
}
\description{
Get reference quality score gamma for each batch
}
