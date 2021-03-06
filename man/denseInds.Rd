% Generated by roxygen2 (4.0.1): do not edit by hand
\name{denseDiagInds}
\alias{denseDiagInds}
\alias{denseLowerInds}
\alias{denseUpperInds}
\title{Indices for square dense matrices}
\usage{
denseDiagInds(n)

denseLowerInds(n)

denseUpperInds(n)
}
\arguments{
\item{n}{matrix size}
}
\description{
Indices for square dense matrices
}
\examples{
set.seed(1)
A <- matrix(rnorm(25), 5, 5)
A[denseDiagInds(5)]  ## = diag(A)
A[denseLowerInds(5)] ## = A[lower.tri(A)]
A[denseUpperInds(5)] ## = A[upper.tri(A)]
}

