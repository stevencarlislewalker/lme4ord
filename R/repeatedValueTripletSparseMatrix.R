##' Sparse matrix with repeated values
##'
##' @param rowInds 1-based row indices (converted to 0-based in
##' output)
##' @param colInds 1-based column indices (converted to 0-based in
##' output)
##' @param valInds indices for the nonzero values
##' @param vals nonzero values
##' @param Dim matrix dimensions
##' @return object of class \code{repSparse} with list elements
##' \code{rowInds}, \code{colInds}, \code{valInds}, \code{vals}, and
##' \code{Dim} attibute
##' @rdname repSparse
##' @export
repSparse <- function(rowInds, colInds, valInds, vals, Dim) {
    if(missing(Dim)) Dim <- c(max(rowInds), max(colInds))
    if(!all(1:max(valInds) %in% valInds)) stop("max(valInds) unnecessarily large")
    if(length(vals) != max(valInds)) stop("mismatch between vals and valInds")
    if(length(rowInds) != length(colInds)) stop("row and column index mismatch")
    if(length(rowInds) != length(valInds)) stop("row and value index mismatch")
    structure(list(rowInds = as.integer(rowInds - 1),
                   colInds = as.integer(colInds - 1),
                   valInds = valInds,
                   vals = vals),
              class = "repSparse",
              Dim = Dim)
}

##' @rdname repSparse
##' @export
print.repSparse <- function(x, ...) str(x, ...)

##' @rdname repSparse
##' @export
as.data.frame.repSparse <- function(x, ...) {
    with(x, {
        data.frame(rowInds = rowInds,
                   colInds = colInds,
                   valInds = valInds,
                   values = vals[valInds])
    })
}

##' @rdname repSparse
##' @param sparse return \code{sparseMatrix}?
##' @export
as.matrix.repSparse <- function(x, sparse = FALSE, ...) {
    ans <- with(x, {
        sparseMatrix(i = rowInds + 1L,
                     j = colInds + 1L,
                     x = vals[valInds],
                     dims = dim(x), ...)
    })
    if(sparse) return(ans)
    return(as.matrix(ans))
}

##' @rdname repSparse
##' @export
dim.repSparse <- function(x) attr(x, "Dim")

##' Convert to \code{sparseMatrix}
##'
##' @param object \code{repSparse} object
##' @return results of \code{sparseMatrix}
##' @rdname repSparse
##' @export
repSparse2Sparse <- function(object, ...) {
    with(object, {
        sparseMatrix(i = rowInds + 1L,
                     j = colInds + 1L,
                     x = vals[valInds],
                     dims = dim(object),
                     giveCsparse = FALSE)
    })
}

##' @rdname repSparse
##' @export
image.repSparse <- function(x, ...) image(repSparse2Sparse(x))


## subset.repSparse <- function(x, rows, cols) {
##     ri <- x$rowInds %in% rows
##     ci <- x$colInds %in% cols
##     inds <- ri & ci
##     print(inds)
##     ans <- within(unclass(x), {
##         rowInds <- rowInds[inds]
##         colInds <- colInds[inds]
##         valInds <- valInds[inds]
##     })
##     with(ans, repSparse(rowInds, colInds, valInds, vals, dim(x)))
## }

## mmult <- function(X, Y) {
    
## }

##' Kronecker and Khatri-Rao products for repeated sparse matrices
##'
##' @param X,Y \code{repSparse} objects
##' @param FUN ignored
##' @param makedimnames ignored
##' @param ... ignored
##' @rdname kron
##' @export
kron <- function(X, Y, FUN = "*",
                 makedimnames = FALSE, ...) {
    lenX <- length(X$rowInds)
    lenY <- length(Y$rowInds)
    lenRep <- rep.int(lenY, lenX)
    XYrowInds <- rep.int(Y$rowInds, lenX) + dim(Y)[1] * rep.int(X$rowInds, lenRep)
    XYcolInds <- rep.int(Y$colInds, lenX) + dim(Y)[2] * rep.int(X$colInds, lenRep)
    XYvalInds <- as.vector(outer(Y$valInds, max(Y$valInds) * (X$valInds - 1), "+"))
    XYvals <- as.vector(Y$vals %*% t(X$vals))
    structure(list(rowInds = XYrowInds,
                   colInds = XYcolInds,
                   valInds = XYvalInds,
                   vals = XYvals),
              class = "repSparse",
              Dim = dim(X) * dim(Y))
}

##' @rdname kron
##' @export
kr <- function(X, Y) {
    matchBin <- outer(X$colInds, Y$colInds, "==")
    matchInd <- which(matchBin, arr.ind = TRUE)

    XrowInds <- X$rowInds[matchInd[, 1]]
    YrowInds <- Y$rowInds[matchInd[, 2]]
    XvalInds <- X$valInds[matchInd[, 1]]
    YvalInds <- Y$valInds[matchInd[, 2]]

    XYcolInds <- X$colInds[matchInd[, 1]]
    XYvalInds <- YvalInds + (length(Y$vals) * (XvalInds - 1))
    XYvals <- as.vector(Y$vals %*% t(X$vals))
    XYrowInds <- YrowInds + (dim(Y)[1] * XrowInds)
    
    structure(list(rowInds = XYrowInds,
                   colInds = XYcolInds,
                   valInds = XYvalInds,
                   vals = XYvals),
              class = "repSparse",
              Dim = c(dim(X)[1] * dim(Y)[1], dim(X)[2]))
}

##' Random repeated sparse matrix
##'
##' @param nrows,ncols numbers of rows and columns
##' @param nvals number of values
##' @param nnonzeros number of nonzero elements
##' @param rfunc random number function
##' @param ... dots
##' @export
rRepSparse <- function(nrows, ncols, nvals, nnonzeros, rfunc = rnorm, ...) {
    if(nnonzeros < nvals)
        stop("number of nonzeros must be at least the number of values")
    if(nnonzeros > (nrows * ncols))
        stop("too many nonzeros for matrix of this size")
    valInds <- sample(c(1:nvals, sample(1:nvals, nnonzeros - nvals, TRUE)))
    indChoose <- sample(nrows * ncols, length(valInds))
    inds <- expand.grid(1:nrows, 1:ncols)[indChoose, ]
    vals <- rfunc(max(valInds), ...)
    repSparse(rowInds = inds$Var1,
              colInds = inds$Var2,
              valInds = valInds, vals = vals,
              Dim = c(nrows, ncols))
}

