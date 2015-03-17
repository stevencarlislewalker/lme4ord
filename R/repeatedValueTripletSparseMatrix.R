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
sparse2RepSparse <- function(object, ...) {
    object <- as(object, "TsparseMatrix")
    repSparse(rowInds = object@i + 1L,
              colInds = object@j + 1L,
              valInds = 1:length(object@i),
              vals = object@x,
              Dim = dim(object))
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

## ----------------------------------------------------------------------
## Matrix multiplication -- mmult (standard matrix product),
## kron (Kronecker product), kr (Khatri-Rao product)
## ----------------------------------------------------------------------

mmult <- function(X, Y) {
    stop("not done")
    matchBin <- outer(X$colInds, Y$rowInds, "==")
    matchInd <- which(matchBin, arr.ind = TRUE)
}

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

##' Simplify a repeated sparse matrix
##'
##' Remove unused values and renumber value indices of a
##' \code{repSparse} object that has values that are not used in any
##' matrix element.
##'
##' @param object a \code{repSparse} object
##' @export
simplifyRepSparse <- function(object) {
    keepers <- sort(unique(object$valInds))
    valsOut <- object$vals[keepers]
    valIndsOut <- match(object$valInds, keepers)
    object$vals <- valsOut
    object$valInds <- valIndsOut
    return(object)
}

## ----------------------------------------------------------------------
## Matrix binding and repeating
## ----------------------------------------------------------------------

##' Row, column, and block-diagonal binding for repeated sparse
##' matrices
##'
##' @param ... list of \code{repSparse} objects
##' @param type type of binding
##' @rdname bind
##' @export
bind <- function(...,
                 type = c("row", "col", "diag")) {
    mats <- list(...)
    nmat <- length(mats)
    if(nmat == 1L) return(mats[[1]])
    with(listTranspose(mats), {
        rowOff <- colOff <- 0
        vilen <- sapply(valInds, length)
        valen <- sapply(vals[-nmat], length) # don't need value lengths for last matrix
        valOff <- rep.int(c(0, cumsum(valen)), vilen)
        if((type == "row") || (type == "diag")) {
            rlen <- sapply(rowInds, length)
            nrows <- sapply(mats, nrow)
            rowOff <- rep.int(c(0, cumsum(nrows[-nmat])), rlen)
        }
        if((type == "col") || (type == "diag")) {
            clen <- sapply(colInds, length)
            ncols <- sapply(mats, ncol)
            colOff <- rep.int(c(0, cumsum(ncols[-nmat])), clen)
        }
        if((type == "row") && (type != "diag")) ncols <- ncol(mats[[1]])
        if((type == "col") && (type != "diag")) nrows <- nrow(mats[[1]])
        ans <- repSparse(rowInds = unlist(rowInds) + rowOff + 1,
                         colInds = unlist(colInds) + colOff + 1,
                         valInds = unlist(valInds) + valOff,
                         vals = unlist(vals),
                         Dim = c(sum(nrows), sum(ncols)))
        return(ans)
    })
}

##' @param lst list of \code{repSparse} objects
##' @rdname bind
##' @export
.bind <- function(lst, type = c("row", "col", "diag")) {
    lapply(bind, c(lst, list(type = type)))
}


##' Repeat repeated sparse matrix
##' @param x \code{repSparse} object
##' @param times like \code{rep}
##' @rdname bind
##' @export
rep.repSparse <- function(x, times,
                          type = c("row", "col", "diag")) {

    len <- length(x$rowInds)
    rowInds <- rep.int(x$rowInds, times)
    colInds <- rep.int(x$colInds, times)
    rowMult <- colMult <- 1
    if((type == "row") || (type == "diag")) {
        off <- rep.int((0:(times-1)) * nrow(x), rep.int(len, times))
        rowInds <- rowInds + off
        rowMult <- times
    }
    if((type == "col") || (type == "diag")) {
        off <- rep.int((0:(times-1)) * ncol(x), rep.int(len, times))
        colInds <- colInds + off
        colMult <- times
    }
    ans <- repSparse(rowInds = rowInds + 1,
                     colInds = colInds + 1,
                     valInds = rep(x$valInds, times = times),
                     vals = x$vals,
                     Dim = c(rowMult, colMult) * dim(x))
    return(ans)
}

## ----------------------------------------------------------------------
## Matrix reshaping -- t
## ----------------------------------------------------------------------

##' Repeated sparse matrix transpose
##'
##' @param x \code{repSparse} object
##' @export
t.repSparse <- function(x) {
    tx <- x
    tx$rowInds <- x$colInds
    tx$colInds <- x$rowInds
    attr(tx, "Dim") <- rev(dim(x))
    return(tx)
}


## ----------------------------------------------------------------------
## Construct special matrices -- repSparseCompSymm
## ----------------------------------------------------------------------

##' Repeated sparse matrix with compound symmetry
##'
##' @param diagVal value for the diagonal
##' @param offDiagVal value for the off-diagonal
##' @param matSize size of the resulting matrix
##' @export
repSparseCompSymm <- function(diagVal, offDiagVal, matSize) {
    iii <- rep.int(1:(matSize-1), 1:(matSize-1)) + 1
    jjj <- sequence(1:(matSize-1))
    ii <- c(1:m, iii, jjj)
    jj <- c(1:m, jjj, iii)
    vi <- rep.int(1:2, c(matSize, 2 * choose(matSize, 2)))
    va <- setNames(c( diagVal ,  offDiagVal ),
                   c("diagVal", "offDiagVal"))
    ans <- repSparse(ii, jj, vi, va, Dim = c(matSize, matSize))
    class(ans) <- c("repSparseCompSymm", class(ans))
    return(ans)
}

##' Diagonal repeated sparse matrix
##'
##' @param vals vector of values
##' @param valInds vector of value indices
##' @export
repSparseDiag <- function(vals, valInds) {
    if(missing(valInds)) valInds <- seq_along(vals)
    matSize <- length(valInds)
    ans <- repSparse(1:matSize, 1:matSize, valInds, vals)
    class(ans) <- c("repSparseDiag", class(ans))
    return(ans)
}

##' Identity repeated sparse matrix
##'
##' @param matSize size of matrix
##' @export
repSparseIdent <- function(matSize) {
    ans <- repSparseDiag(1, rep(1, matSize))
    class(ans) <- c("repSparseIdent", class(ans))
    return(ans)
}
