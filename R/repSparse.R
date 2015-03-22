##' Sparse matrix with repeated values
##'
##' @param rowInds 1-based row indices (converted to 0-based in
##' output)
##' @param colInds 1-based column indices (converted to 0-based in
##' output)
##' @param valInds indices for the nonzero values
##' @param vals nonzero values
##' @param trans function for transforming some parameters into
##' \code{vals}.  the \code{\link{environment}} of this function must
##' contain a named list of these parameters (called \code{matPars}),
##' which are to be taken as the arguments to \code{trans}.
##' \code{update} will update these parameters (in the future).
##' @param Dim matrix dimensions
##' @return object of class \code{repSparse} with list elements
##' \code{rowInds}, \code{colInds}, \code{valInds}, \code{vals}, and
##' \code{trans}, and a \code{Dim} attibute
##' @rdname repSparse
##' @export
repSparse <- function(rowInds, colInds, valInds, vals, trans, Dim) {
    if(missing(Dim)) Dim <- c(max(rowInds), max(colInds))
    if(!all(1:max(valInds) %in% valInds)) stop("max(valInds) unnecessarily large")
    if(length(vals) != max(valInds)) stop("mismatch between vals and valInds")
    if(length(rowInds) != length(colInds)) stop("row and column index mismatch")
    if(length(rowInds) != length(valInds)) stop("row and value index mismatch")
    if(missing(trans)) trans <- mkIdentityTrans()
    structure(list(rowInds = as.integer(rowInds - 1),
                   colInds = as.integer(colInds - 1),
                   valInds = valInds,
                   vals = vals,
                   trans = trans),
              class = "repSparse",
              Dim = Dim)
}


##' @rdname repSparse
##' @export
print.repSparse <- function(x, ...) str(x, ...)

##' @param newPars new parameter values
##' @rdname repSparse
##' @export
update.repSparse <- function(object, newPars, ...) {
    object$vals <- object$trans(newPars)
    return(object)
}


##' @param decreasing see \code{\link{sort}}
##' @param type sort by column, row, or value indices?
##' @rdname repSparse
##' @export
sort.repSparse <- function(x, decreasing = FALSE,
                           type = c("col", "row", "val"), ...) {
    inds <- switch(type[1],
                   col = x$colInds,
                   row = x$rowInds,
                   val = x$valInds)
    ord <- order(inds, decreasing = decreasing, ...)
    x$rowInds <- x$rowInds[ord]
    x$colInds <- x$colInds[ord]
    x$valInds <- x$valInds[ord]
    return(x)
}


##' @rdname repSparse
##' @export
as.data.frame.repSparse <- function(x, ...) {
    with(x, {
        data.frame(rowInds = rowInds,
                   colInds = colInds,
                   valInds = valInds,
                   vals = vals[valInds])
    })
}

##' @param keep.rownames passed to \code{as.data.table}
##' @rdname repSparse
##' @export
as.data.table.repSparse <- function(x, keep.rownames = FALSE) {
    as.data.table(as.data.frame(x), keep.rownames = keep.rownames)
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
    x <- object@x
    
                                        # try to find detectable
                                        # repeated structure.
    va <- sort(unique(x))
    vi <- match(x, va)

                                        # return value
    repSparse(rowInds = object@i + 1L,
              colInds = object@j + 1L,
              valInds = vi,
              vals = va,
              trans = mkIdentityTrans(),
              Dim = dim(object))
}


##' @rdname repSparse
##' @export
image.repSparse <- function(x, ...) image(repSparse2Sparse(x))

##' @param rows,cols not sure yet
##' @rdname repSparse
##' @export
subset.repSparse <- function(x, rows, cols) {
    stop("not done, but really should be")
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
}



## ----------------------------------------------------------------------
## Matrix operations -- mmult (standard matrix product), kron
## (Kronecker product), kr (Khatri-Rao product), madd (standard matrix
## addition)
## ----------------------------------------------------------------------


##' Row-wise combination of two repeated sparse matrices
##'
##' @param X,Y \code{repSparse} objects
##' @return a row-wise combination of repeated sparse matrices
##' @family matrixOperations
##' @export
rowWiseCombination <- function(X, Y) {

    matchBin <- outer(X$colInds, Y$colInds, "==") ## this is slow!
    matchInd <- which(matchBin, arr.ind = TRUE)
    
    structure(list(XrowInds  = X$rowInds[matchInd[, 1]],
                   YrowInds  = Y$rowInds[matchInd[, 2]],
                   XYcolInds = X$colInds[matchInd[, 1]],
                   XvalInds  = X$valInds[matchInd[, 1]],
                   YvalInds  = Y$valInds[matchInd[, 2]],
                   Xvals = X$vals,
                   Yvals = Y$vals,
                   Xtrans = X$trans,
                   Ytrans = Y$trans),
              class = "repSparseRowWiseCombination",
              Dim = c(nrow(X), nrow(Y), ncol(X)))
}



##' Standard matrix multiplication for repeated sparse matrices
##'
##' @param X,Y \code{repSparse} objects
##' @param trans two argument transformation from \code{X$vals} and
##' \code{Y$vals} to the repeated values of the resulting matrix
##' @family matrixOperations
##' @export
mmult <- function(X, Y, trans = "*") {

    with(rowWiseCombination(X, t(Y)), {

        outerColInds <- XYcolInds
        outerValInds <- YvalInds + (length(Yvals) * (XvalInds - 1))
        outerTrans <- mkOuterTrans(Xvals, Yvals, Xtrans, Ytrans, trans)
        outerVals <- as.vector(outer(Xvals, Yvals, trans))
        outerRowInds <- YrowInds + (dim(Y)[2] * XrowInds)
        
        sumInds <- lapply(tapply(outerValInds, outerRowInds, I), sort)
        uniqueSumInds <- unique(sumInds)
        
        groupNames <- as.numeric(names(sumInds))
        matchNames <- match(groupNames, outerRowInds)
        
        XYrowInds <- XrowInds[matchNames]
        XYcolInds <- YrowInds[matchNames]
        XYvals <- sapply(uniqueSumInds, function(ii) sum(outerVals[ii]))
        XYvalInds <- match(sumInds, uniqueSumInds)
        
        structure(list(rowInds = XYrowInds,
                       colInds = XYcolInds,
                       valInds = XYvalInds,
                       vals = XYvals,
                       trans = outerTrans),
                  class = "repSparse",
                  Dim = c(dim(X)[1], dim(Y)[2]))
    })
}


##' Kronecker and Khatri-Rao products for repeated sparse matrices
##'
##' @param X,Y \code{repSparse} objects
##' @param trans see argument \code{FUN} in \code{\link{outer}}
##' @param makedimnames ignored
##' @param ... ignored
##' @rdname kron
##' @family matrixOperations
##' @export
kron <- function(X, Y, trans = "*",
                 makedimnames = FALSE, ...) {
    lenX <- length(X$rowInds)
    lenY <- length(Y$rowInds)
    lenRep <- rep.int(lenY, lenX)
    XYrowInds <- rep.int(Y$rowInds, lenX) + dim(Y)[1] * rep.int(X$rowInds, lenRep)
    XYcolInds <- rep.int(Y$colInds, lenX) + dim(Y)[2] * rep.int(X$colInds, lenRep)
    XYvalInds <- as.vector(outer(Y$valInds, max(Y$valInds) * (X$valInds - 1), "+"))
    XYtrans <- mkOuterTrans(Y$vals, X$vals, Y$trans, X$trans, trans)
    XYvals <- as.vector(outer(Y$vals, X$vals, FUN = trans))
    structure(list(rowInds = XYrowInds,
                   colInds = XYcolInds,
                   valInds = XYvalInds,
                   vals = XYvals,
                   trans = XYtrans),
              class = "repSparse",
              Dim = dim(X) * dim(Y))
}

##' @rdname kron
##' @family matrixOperations
##' @export
kr <- function(X, Y, trans = "*") {

    with(rowWiseCombination(X, Y), {
    
        XYvalInds <- YvalInds + (length(Yvals) * (XvalInds - 1))
        XYtrans <- mkOuterTrans(Yvals, Xvals, Ytrans, Xtrans, trans)
        XYvals <- as.vector(outer(Yvals, Xvals, FUN = trans))
        XYrowInds <- YrowInds + (dim(Y)[1] * XrowInds)
        
        structure(list(rowInds = XYrowInds,
                       colInds = XYcolInds,
                       valInds = XYvalInds,
                       vals = XYvals,
                       trans = XYtrans),
                  class = "repSparse",
                  Dim = c(dim(X)[1] * dim(Y)[1], dim(X)[2]))
    })
}



##' Simplify a repeated sparse matrix
##'
##' Remove unused values and renumber value indices of a
##' \code{repSparse} object that has values that are not used in any
##' matrix element.  (FIXME: not currently used anywhere)
##'
##' @param object a \code{repSparse} object
##' @export
simplifyRepSparse <- function(object) {
    warning("untested")
    keepers <- sort(unique(object$valInds))
    valsOut <- object$vals[keepers]
    valIndsOut <- match(object$valInds, keepers)
    object$vals <- valsOut
    object$valInds <- valIndsOut
    return(object)
}



## ----------------------------------------------------------------------
## Make trans functions
## ----------------------------------------------------------------------

##' Construct functions for transforming a parameter vector to the
##' non-zero values of a repeated sparse matrix
##'
##' @rdname mkTrans
##' @family mkTransFunctions
##' @export
mkIdentityTrans <- function() {
    return(function(matPars) matPars)
}


##' @param A,B \code{numeric} \code{vector}s
##' @param Atrans,Btrans functions for transforming \code{A} and
##' \code{B}
##' @param ABtrans function to pass as \code{FUN} in
##' \code{\link{outer}}
##' @rdname mkTrans
##' @family mkTransFunctions
##' @export
mkOuterTrans <- function(A, B, Atrans, Btrans, ABtrans) {

                                        # check to see if one of the
                                        # parameter sets is a single
                                        # 1L, in which case just
                                        # return an identity
                                        # transformation (because an
                                        # outer product involving a
                                        # scalar of 1L is boring)
    checkForBordom <- function(xx) (length(xx) == 1L) & (xx[1] == 1L)
    if(checkForBordom(A) || checkForBordom(B)) return(mkIdentityTrans())

                                        # construct the environment of
                                        # the transformation function
    local({
        ## FIXME: include indices??  don't think so, but ... actually
        ## now i think yes
        Atrans <- Atrans
        Btrans <- Btrans
        ABtrans <- ABtrans
        Aind <- seq_along(A)
        Bind <- seq_along(B) + length(A)
        function(matPars) as.vector(outer(Atrans(matPars[Aind]),
                                          Btrans(matPars[Bind]),
                                          FUN = ABtrans))
    })
}

##' @param valsList list of value vectors
##' @param transList list of transform functions
##' @rdname mkTrans
##' @family mkTransFunctions
##' @export
mkListTrans <- function(valsList, transList) {
    local({
        transList <- transList
        indList <- lapply(valsList, seq_along)
        function(matPars) {
            unlist(lapply(seq_along(indList),
                          function(ii) transList[[i]](matPars[indList[[i]]])))
        }
    })
}

##' @rdname mkTrans
##' @family mkTransFunctions
##' @export
mkCholOneOffDiagTrans <- function() {
    function(matPars) {
        with(setNames(as.list(matPars), c("diagVal", "offDiagVal")), {
            c(sqrt(diagVal),
              offDiagVal/sqrt(diagVal),
              sqrt(diagVal - ((offDiagVal^2) / diagVal)))
        })
    }
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
##' @family matrixBinding
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
                         trans = mkListTrans(vals, trans),
                         Dim = c(sum(nrows), sum(ncols)))
        return(ans)
    })
}

##' @param lst list of \code{repSparse} objects
##' @rdname bind
##' @family matrixBinding
##' @export
.bind <- function(lst, type = c("row", "col", "diag")) {
    lapply(bind, c(lst, list(type = type)))
}


##' Repeat repeated sparse matrix
##' @param x \code{repSparse} object
##' @param times like \code{rep}
##' @rdname bind
##' @family matrixBinding
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
                     trans = x$trans,
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
## Changing sparse formats 
## ----------------------------------------------------------------------

##' Changing sparse format
##'
##' @param point vector of column pointers
##' @rdname changeSparseFormat
##' @export
point2ind <- function(point) {
                                        # ?Matrix::sparseMatrix
    dp <- diff(sort(point))
    rep.int(seq_along(dp), dp)
}

##' @param ind vector of column indices
##' @param maxInd number of rows or columns (could be larger than
##' \code{max(ind)})
##' @param fillNA fill in \code{NA}'s with repeated column numbers (as
##' in \code{Matrix} package)?
##' @rdname changeSparseFormat
##' @export
ind2point <- function(ind, maxInd, fillNA = TRUE) {
    point <- match(0:maxInd, ind) - 1L
    if(!fillNA) return(point)
    point[maxInd + 1] <- maxInd
    isNaPoint <- is.na(point)
    for(j in maxInd:1) {
        point[j] <- ifelse(isNaPoint[j],
                           point[j + 1],
                           point[j])
    }
    return(point)
}

## ----------------------------------------------------------------------
## Construct special matrices -- repSparseCompSymm, repSparseDiag,
## repSparseIdent, rRepSparse
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
    ii <- c(1:m, iii, jjj) ## FIXME: where did m come from?!
    jj <- c(1:m, jjj, iii)
    vi <- rep.int(1:2, c(matSize, 2 * choose(matSize, 2)))
    va <- setNames(c( diagVal ,  offDiagVal ),
                   c("diagVal", "offDiagVal"))
    ans <- repSparse(ii, jj, vi, va, Dim = c(matSize, matSize))
    class(ans) <- c("repSparseCompSymm", class(ans))
    return(ans)
}

##' Repeated sparse covariance matrix with equal variances and one
##' off-diagonal covariance
##'
##' @param diagVal value for the diagonal
##' @param offDiagVal value for the off-diagonal
##' @param offDiagInds indices for the two correlated objects
##' @param matSize size of the resulting matrix
##' @export
repSparseOneOffDiag <- function(diagVal, offDiagVal, offDiagInds, matSize) {
    if(length(offDiagInds) != 2L) stop("only one off diagonal element please")
    if(offDiagInds[1] == offDiagInds[2]) stop("off diagonal must be off the diagonal")
    iii <- jjj <- 1:matSize
    ii <- c(iii,     offDiagInds )
    jj <- c(jjj, rev(offDiagInds))
    vi <- c(rep(1, matSize), rep(2, 2))
    va <- setNames(c( diagVal ,  offDiagVal ),
                   c("diagVal", "offDiagVal"))
    ans <- repSparse(ii, jj, vi, va, Dim = c(matSize, matSize))
    class(ans) <- c("repSparseOneOffDiag", class(ans))
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


## ----------------------------------------------------------------------
## Cholesky -- 
## ----------------------------------------------------------------------

##' Cholesky decomposition of repeated sparse matrices
##'
##' @note These are often just bailout methods, but some special
##' \code{repSparse} matrices have exploitable structure, which can be
##' used to keep the number of repeated values down.
##' @param x an object that inherits from class
##' \code{\link{repSparse}}
##' @param ... passed to subsequent functions
##' @rdname chol
##' @export
chol.repSparse <- function(x, ...) {
    sparse2RepSparse(chol(as.matrix(x, sparse = TRUE)))
}

##' @rdname chol
##' @export
chol.repSparseOneOffDiag <- function(x, ...) {
    offRow <- sort(x$rowInds[c(0, -1) + length(x$valInds)]) + 1L
    va <- c(sqrt(x$vals[1]),
            x$vals[2]/sqrt(x$vals[1]),
            sqrt(x$vals[1] - ((x$vals[2]^2) / x$vals[1])))
    ni <- length(x$valInds)
    vi <- x$valInds[-ni]
    ri <- x$rowInds[-ni] + 1L
    ci <- x$colInds[-ni] + 1L
    vi[offRow[2]] <- 3
    vi[length(vi)] <- 2
    ans <- repSparse(ri, ci, vi, va,
                     trans = mkCholOneOffDiagTrans(),
                     Dim = dim(x))
    return(ans)
}

