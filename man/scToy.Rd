% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SoupX.R
\docType{data}
\name{scToy}
\alias{scToy}
\title{Toy SoupChanel object}
\format{
\code{scToy} is a \code{SoupChannel} object.
}
\usage{
data(scToy)
}
\description{
A \code{\link{SoupChannel}} object created from the toy data used in examples.
}
\details{
The toy data is created from a modified version of the extremely reduced \code{Seurat} \code{pbmc_small} dataset.  It includes clusters, tSNE coordinates and a flat estimate of 0.1 contamination.  It includes data for only 226 genes and 62 cells and should not be used for anything other than testing functions as it is not representative of real data in any way.
}
\keyword{datasets}
