##' Sparse matrix with repeated values
##'
##' @param rowInds 1-based row indices (converted to 0-based in
##' output)
##' @param colInds 1-based column indices (converted to 0-based in
##' output)
##' @param valInds indices for the nonzero values
##' @param vals nonzero values
##' @param trans function for transforming some parameters into
##' \code{vals}.  the \code{\link{environment}} of this function
##' should contain a vector of these parameters (called \code{init}),
##' which could be taken as the arguments to \code{trans}.  these
##' initial values can be obtained by \code{\link{getInit}} and set by
##' \code{\link{setInit}}.  an \code{\link{update.repSparse}} method
##' will update these parameters.
##' @param Dim matrix dimensions
##' @return object of class \code{repSparse} with list elements
##' \code{rowInds}, \code{colInds}, \code{valInds}, \code{vals}, and
##' \code{trans}, and a \code{Dim} attibute
##' @rdname repSparse
##' @family repSparseTopics
##' @export
repSparse <- function(rowInds, colInds, valInds, vals, trans, Dim) {
    if(missing(Dim)) Dim <- c(max(rowInds), max(colInds))
    if(!all(1:max(valInds) %in% valInds))  stop("max(valInds) unnecessarily large")
    if(length(vals) != max(valInds))       stop("mismatch between vals and valInds")
    if(length(rowInds) != length(colInds)) stop("row and column index mismatch")
    if(length(rowInds) != length(valInds)) stop("row and value index mismatch")
    if(missing(trans)) trans <- mkIdentityTrans(vals)
    structure(list(rowInds = as.integer(rowInds - 1),
                   colInds = as.integer(colInds - 1),
                   valInds = valInds,
                   vals = vals,
                   trans = trans),
              class = "repSparse",
              Dim = Dim)
}


## ----------------------------------------------------------------------
## repSparse-class
## ----------------------------------------------------------------------

##' \code{repSparse} class
##'
##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparse
setOldClass("repSparse")

##' @param x \code{repSparse} object
##' @param ... passed to subsequent functions
##' @rdname repSparse-class
##' @export
print.repSparse <- function(x, ...) {
    attr(x, "components") <- NULL
    str(x, ...)
}

##' @param object \code{repSparse} object
##' @param newPars new parameter values
##' @rdname repSparse-class
##' @export
update.repSparse <- function(object, newPars, ...) {
    if(missing(newPars)) newPars <- getInit(object)
    object$vals <- object$trans(newPars)
    return(object)
}

##' @rdname repSparse-class
##' @export
image.repSparse <- function(x, ...) image(as.matrix(x, sparse = TRUE))

##' @rdname repSparse-class
##' @export
t.repSparse <- function(x) {
    tx <- x
    tx$rowInds <- x$colInds
    tx$colInds <- x$rowInds
    attr(tx, "Dim") <- rev(dim(x))
    return(tx)
}

##' @param decreasing see \code{\link{sort}}
##' @param type sort by column, row, or value indices?
##' @rdname repSparse-class
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

##' @rdname repSparse-class
##' @export
dim.repSparse <- function(x) attr(x, "Dim")


## ----------------------------------------------------------------------
## Coercion -- as...
## ----------------------------------------------------------------------

##' Coerce to and from repeated sparse matrices
##'
##' @param x an object
##' @param ... dots
##' @rdname as.repSparse
##' @family repSparseTopics
##' @export
##' @examples
##' set.seed(1)
##' m1 <- as.repSparse(matrix(rnorm(6), 2, 3))
##' m2 <- as.repSparse(matrix(rbinom(6, 1, 0.5), 2, 3))
##' m3 <- sparseMatrix(i = 1:10, j = rep(1:5, 2),
##'                    x = as.integer(rbinom(10, 1, 0.5)),
##'                    giveCsparse = FALSE)
##' as.repSparse(m1)
##' as.repSparse(m2)
##' as.repSparse(m3)
##' as(m1, "TsparseMatrix")
##' as(m2, "dgCMatrix")
##' as(m2, "sparseMatrix")
as.repSparse <- function(x, ...) {
    UseMethod("as.repSparse")
}

##' @rdname as.repSparse
##' @export
as.repSparse.repSparse <- function(x, ...) {
    class(x) <- "repSparse"
    return(x)
}

##' @rdname as.repSparse
##' @export
as.repSparse.dsparseMatrix <- function(x, ...) {
    x <- as(x, "TsparseMatrix")
                                        # try to find detectable
                                        # repeated structure.
    va <- sort(unique(x@x))
    vi <- match(x@x, va)

                                        # return value
    ans <- repSparse(rowInds = x@i + 1L,
                     colInds = x@j + 1L,
                     valInds = vi,
                     vals = va,
                     trans = mkIdentityTrans(va),
                     Dim = dim(x))
    return(sort(ans))
}

##' @rdname as.repSparse
##' @export
as.repSparse.matrix <- function(x, ...) {
    vecx <- as.numeric(x)
    va <- sort(unique(vecx))
    vi <- match(vecx, va)
    ii <- rep.int(1:nrow(x), ncol(x))
    jj <- rep.int(1:ncol(x), rep.int(nrow(x), ncol(x)))
    ans <- repSparse(rowInds = ii, colInds = jj,
                     valInds = vi, vals = va,
                     trans = mkIdentityTrans(va),
                     Dim = dim(x))
    return(sort(ans))
}

##' @rdname as.repSparse
##' @export
as.repSparse.factor <- function(x, ...) {
    as.repSparse(as(x, "sparseMatrix"))
}

##' @rdname as.repSparse
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
##' @rdname as.repSparse
##' @export
as.data.table.repSparse <- function(x, keep.rownames = FALSE) {
    as.data.table(as.data.frame(x), keep.rownames = keep.rownames)
}

##' @rdname as.repSparse
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

repSparse2Csparse <- function(from) {
    ans <- new("dgCMatrix")
    ans@i <- as.integer(from$rowInds)
    ans@p <- as.integer(ind2point(from$colInds, ncol(from)))
    ans@x <- as.numeric(with(from, vals[valInds]))
    ans@Dim <- dim(from)
    return(ans)
}
repSparse2Tsparse <- function(from) {
    ans <- new("dgTMatrix")
    ans@i <- as.integer(from$rowInds)
    ans@j <- as.integer(from$colInds)
    ans@x <- as.numeric(with(from, vals[valInds]))
    ans@Dim <- dim(from)
    return(ans)
}

##' as("repSparse", "sparseMatrix")
##' @name as
##' @rdname as.repSparse
##' @importClassesFrom Matrix sparseMatrix
setAs("repSparse",  "sparseMatrix", def = repSparse2Csparse)

##' as("repSparse", "CsparseMatrix")
##' @name as
##' @rdname as.repSparse
##' @importClassesFrom Matrix CsparseMatrix
setAs("repSparse", "CsparseMatrix", def = repSparse2Csparse)

##' as("repSparse", "dgCMatrix")
##' @name as
##' @rdname as.repSparse
##' @importClassesFrom Matrix dgCMatrix
setAs("repSparse",     "dgCMatrix", def = repSparse2Csparse)

##' as("repSparse", "TsparseMatrix")
##' @name as
##' @rdname as.repSparse
##' @importClassesFrom Matrix TsparseMatrix
setAs("repSparse", "TsparseMatrix", def = repSparse2Tsparse)

##' as("repSparse", "dgTMatrix")
##' @name as
##' @rdname as.repSparse
##' @importClassesFrom Matrix dgTMatrix
setAs("repSparse",     "dgTMatrix", def = repSparse2Tsparse)


## ----------------------------------------------------------------------
## Initial values -- get and set init parameter vectors for repeated
## sparse matrices
## ----------------------------------------------------------------------

##' Get and set initial parameter values for repeated sparse matrices
##'
##' @param x object
##' @param ... not yet used
##' @rdname getInit
##' @family repSparseTopics
##' @export
##' @examples
##' set.seed(1)
##' m1 <- as.repSparse(matrix(rnorm(6), 2, 3))
##' m2 <- repSparseCompSymm(1.2, -0.2, 5)
##' getInit(m1)
##' getInit(m2)
getInit <- function(x, ...) UseMethod("getInit")

##' @rdname getInit
##' @export
getInit.default <- function(x, ...) x$init

##' @rdname getInit
##' @export
getInit.repSparse <- function(x, ...) environment(x$trans)$init

##' @rdname getInit
##' @export
getInit.function <- function(x, ...) environment(x)$init

##' @param init initial value for the parameter vector
##' @rdname getInit
##' @export
setInit <- function(x, init, ...) UseMethod("setInit")

##' @rdname getInit
##' @export
setInit.default <- function(x, init, ...) assign("init", init, envir = as.environment(x))

##' @rdname getInit
##' @export
setInit.repSparse <- function(x, init, ...) assign("init", init, envir = environment(x$trans))

##' @rdname getInit
##' @export
setInit.function <- function(x, init, ...) assign("init", init, envir = environment(x))


## ----------------------------------------------------------------------
## Matrix operations -- mmult (standard matrix product), kron
## (Kronecker product), kr (Khatri-Rao product), madd (standard matrix
## addition)
## ----------------------------------------------------------------------

##' Row-wise combination of two repeated sparse matrices
##'
##' @param X,Y \code{repSparse} objects
##' @return a row-wise combination of repeated sparse matrices
##' @rdname matrixOperations
##' @family repSparseTopics
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
##' @rdname matrixOperations
##' @export
mmult <- function(X, Y, trans = "*") {

    warning("this is slow and bad and probably even just wrong.  ",
            "please construct standard matrix products with standard tools.")

    with(rowWiseCombination(X, t(Y)), {

        outerColInds <- XYcolInds
        outerValInds <- YvalInds + (length(Yvals) * (XvalInds - 1))
        outerTrans <- mkOuterTrans(Xtrans, Ytrans, trans)
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
##' @param trans see argument \code{FUN} in \code{\link{outer}}
##' @param makedimnames ignored
##' @param saveComponents should component matrices be saved?
##' @param ... ignored
##' @rdname matrixOperations
##' @export
kron <- function(X, Y, trans = "*",
                 makedimnames = FALSE, saveComponents = FALSE, ...) {

    if(saveComponents) {
        components <- list(FUN = kron, X = X, Y = Y, trans = trans)
    } else {
        components <- NULL
    }
    
    lenX <- length(X$rowInds)
    lenY <- length(Y$rowInds)
    lenRep <- rep.int(lenY, lenX)
    XYrowInds <- rep.int(Y$rowInds, lenX) + dim(Y)[1] * rep.int(X$rowInds, lenRep)
    XYcolInds <- rep.int(Y$colInds, lenX) + dim(Y)[2] * rep.int(X$colInds, lenRep)
    XYvalInds <- as.vector(outer(Y$valInds, max(Y$valInds) * (X$valInds - 1), "+"))
    XYtrans <- mkOuterTrans(Y$trans, X$trans, trans)
    XYvals <- as.vector(outer(Y$vals, X$vals, FUN = trans))
    structure(list(rowInds = XYrowInds,
                   colInds = XYcolInds,
                   valInds = XYvalInds,
                   vals = XYvals,
                   trans = XYtrans),
              class = c("repSparseKron", "repSparse"),
              Dim = dim(X) * dim(Y),
              components = components)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseKron
setOldClass("repSparseKron")
setIs("repSparseKron", "repSparse")

##' @rdname matrixOperations
##' @export
kr <- function(X, Y, trans = "*", saveComponents = FALSE) {

    if(saveComponents) {
        components <- list(FUN = kr, X = X, Y = Y, trans = trans)
    } else {
        components <- NULL
    }
    
    with(rowWiseCombination(X, Y), {
        
        XYvalInds <- YvalInds + (length(Yvals) * (XvalInds - 1))
        XYtrans <- mkOuterTrans(Ytrans, Xtrans, trans)
        XYvals <- as.vector(outer(Yvals, Xvals, FUN = trans))
        XYrowInds <- YrowInds + (dim(Y)[1] * XrowInds)
        
        structure(list(rowInds = XYrowInds,
                       colInds = XYcolInds,
                       valInds = XYvalInds,
                       vals = XYvals,
                       trans = XYtrans),
                  class = c("repSparseKr", "repSparse"),
                  Dim = c(dim(X)[1] * dim(Y)[1], dim(X)[2]),
                  components = components)
    })
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseKr
setOldClass("repSparseKr")
setIs("repSparseKr", "repSparse")


## ----------------------------------------------------------------------
## Make trans functions
## ----------------------------------------------------------------------

##' Construct functions for transforming a parameter vector to the
##' non-zero values of a repeated sparse matrix
##'
##' Each \code{trans} function must take one vector argument called
##' \code{matPars}, which are the parameters of the matrix.  The
##' environment of these transformation functions must contain a
##' vector called \code{init}, which contains the initial values (aka,
##' the prototype) of \code{matPars}.
##'
##' @rdname mkTrans
##' @family repSparseTopics
##' @export
mkIdentityTrans <- function(init) {
    local({
        init <- init
        function(matPars) matPars
    })
}

##' @param Atrans,Btrans functions for transforming \code{A} and
##' \code{B}
##' @param ABtrans function to pass as \code{FUN} in
##' \code{\link{outer}}
##' @rdname mkTrans
##' @export
mkOuterTrans <- function(Atrans, Btrans, ABtrans) {
    A <- Atrans(getInit(Atrans))
    B <- Btrans(getInit(Btrans))

                                        # check to see if one of the
                                        # parameter sets is a single
                                        # 1L, in which case just
                                        # return an identity
                                        # transformation (because an
                                        # outer product involving a
                                        # scalar of 1L is boring)
    checkForBordom <- function(xx) (length(xx) == 1L) & (xx[1] == 1L)
    ## if(checkForBordom(A) && !checkForBordom(B)) return(mkIdentityTrans(getInit(Btrans)))
    ## if(!checkForBordom(A) && checkForBordom(B)) return(mkIdentityTrans(getInit(Atrans)))
    ## if(checkForBordom(A) && checkForBordom(B)) return(mkIdentityTrans(1)) ## FIXME: too
                                                                             ## presumptuous?

                                        # construct the environment of
                                        # the transformation function
    local({
        Atrans <- Atrans
        Btrans <- Btrans
        ABtrans <- ABtrans
        Ainit <- getInit(Atrans)
        Binit <- getInit(Btrans)
        Aind <- seq_along(Ainit)
        Bind <- seq_along(Binit) + length(Ainit)
        init <- c(Ainit, Binit)
        function(matPars) as.vector(outer(Atrans(matPars[Aind]),
                                          Btrans(matPars[Bind]),
                                          FUN = ABtrans))
    })
}

##' @param transList list of transform functions
##' @rdname mkTrans
##' @export
mkListTrans <- function(transList) {
    local({
        transList <- transList
        parList <- lapply(transList, getInit)
        indList <- lapply(parList, seq_along)
        for(ii in 2:length(indList)) {
            indList[[ii]] <- indList[[ii]] + max(0L, indList[[ii-1]])
        }
        init <- unlist(parList)
        function(matPars) {
            unlist(lapply(seq_along(indList),
                          function(ii) transList[[ii]](matPars[indList[[ii]]])))
        }
    })
}

##' @param init initial parameter values
##' @rdname mkTrans
##' @export
mkCholOneOffDiagTrans <- function(init) {
    local({
        init <- init
        function(matPars) {
            with(setNames(as.list(matPars), c("diagVal", "offDiagVal")), {
                c(sqrt(diagVal),
                  offDiagVal/sqrt(diagVal),
                  sqrt(diagVal - ((offDiagVal^2) / diagVal)))
            })
        }
    })
}


##' @rdname mkTrans
##' @export
mkGenCholInPlaceTrans <- function(init) {
    local({
        init <- init
        m <- length(init)
        n <- as.integer((sqrt(1 + 8 * m)-1)/2)
        diagInds <- choose(3:(n + 1), 2)
        col1Inds <- choose(1:n, 2) + 1L
        function(matPars) {
            print(matPars[col1Inds] <- matPars[col1Inds]/sqrt(matPars[1]))
            matPars
            k <- 3
            for(i in 3:n) {
                for(j in 1:(i-1)) {
                    k <- k + 1
                    matPars[k] <-
                        (matPars[k] - sum(matPars[col1Inds[i-1] + (0:(j-1))] *
                                          matPars[col1Inds[j-1] + (0:(j-1))])) /
                        matPars[diagInds[i]]
                }
                k <- k + 1
                matPars[k] <- sqrt(matPars[k] - sum(matPars[(k-i-1):(k-1)]^2))
            }
            return(matPars)
        }
    })
}


## ----------------------------------------------------------------------
## Reset trans functions
## ----------------------------------------------------------------------

##' Reset the transformation function of a repeated sparse matrix
##'
##' @param object repeated sparse matrix
##' @rdname resetTrans
##' @export
resetTransConst <- function(object) {
    object$trans <- local({
        baseline <- object$vals
        init <- numeric(0)
        function(matPars = NULL) {
            return(baseline)
        }
    })
    return(object)
}



## ----------------------------------------------------------------------
## Matrix binding and repeating
## ----------------------------------------------------------------------

##' Row, column, and block-diagonal binding for repeated sparse
##' matrices
##'
##' @param ... list of \code{repSparse} objects (but not used for
##' \code{rep.repSparse})
##' @param type type of binding
##' @param saveComponents should component matrices be saved?
##' @rdname bind
##' @family repSparseTopics
##' @export
bind <- function(...,
                 type = c("row", "col", "diag"),
                 saveComponents = FALSE) {
    mats <- list(...)

    if(saveComponents) {
        components <- list(FUN = .bind, mats = mats, type = type)
    } else {
        components <- NULL
    }
    
    nmat <- length(mats)
    type <- type[[1]]
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
                         trans = mkListTrans(trans),
                         Dim = c(sum(nrows), sum(ncols)))
        class(ans) <- c("repSparseBind", class(ans))
        ans$components <- components
        return(ans)
    })
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseBind
setOldClass("repSparseBind")
setIs("repSparseBind", "repSparse")

##' @param mats list of \code{repSparse} matrix objects
##' @rdname bind
##' @export
.bind <- function(mats, type = c("row", "col", "diag")) {
    do.call(bind, c(mats, list(type = type)))
}

##' Repeat repeated sparse matrix
##' @param x \code{repSparse} object
##' @param times like \code{rep}
##' @rdname bind
##' @export
rep.repSparse <- function(x, times,
                          type = c("row", "col", "diag"),
                          ...) {

    type = type[[1]]
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
    class(ans) <- c("repSparseRep", class(ans))
    ans$components <- list(FUN = rep, x = x, times = times, type = type)
    return(ans)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseRep
setOldClass("repSparseRep")
setIs("repSparseRep", "repSparse")


## ----------------------------------------------------------------------
## Changing sparse formats 
## ----------------------------------------------------------------------

##' Changing sparse format
##'
##' @param point vector of column pointers
##' @rdname changeSparseFormat
##' @family repSparseTopics
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
    point[maxInd + 1] <- max(maxInd, length(ind))
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
##' @rdname specialRepSparse
##' @family repSparseTopics
##' @export
repSparseCompSymm <- function(diagVal, offDiagVal, matSize) {
    if((!(diagVal > offDiagVal)) || (!(offDiagVal > (-diagVal)/(matSize-1))))
        warning("resulting matrix not positive definite")
    iii <- rep.int(1:(matSize-1), 1:(matSize-1)) + 1L
    jjj <- sequence(1:(matSize-1))
    ii <- c(1:matSize, iii, jjj)
    jj <- c(1:matSize, jjj, iii)
    vi <- rep.int(1:2, c(matSize, 2 * choose(matSize, 2)))
    va <- setNames(c( diagVal ,  offDiagVal ),
                   c("diagVal", "offDiagVal"))
    ans <- repSparse(ii, jj, vi, va, Dim = c(matSize, matSize))
    class(ans) <- c("repSparseCompSymm", class(ans))
    return(ans)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseCompSymm
setOldClass("repSparseCompSymm")
setIs("repSparseCompSymm", "repSparse")

##' @param offDiagInds indices for the two correlated objects
##' @rdname specialRepSparse
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

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseOneOffDiag
setOldClass("repSparseOneOffDiag")
setIs("repSparseOneOffDiag", "repSparse")

##' Diagonal repeated sparse matrix
##'
##' @param vals vector of values
##' @param valInds vector of value indices
##' @rdname specialRepSparse
##' @export
repSparseDiag <- function(vals, valInds) {
    if(missing(valInds)) valInds <- seq_along(vals)
    matSize <- length(valInds)
    ans <- repSparse(1:matSize, 1:matSize, valInds, vals)
    class(ans) <- c("repSparseDiag", class(ans))
    return(ans)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseDiag
setOldClass("repSparseDiag")
setIs("repSparseDiag", "repSparse")

##' Identity repeated sparse matrix
##'
##' @rdname specialRepSparse
##' @export
repSparseIdent <- function(matSize) {
    ans <- repSparseDiag(1, rep(1, matSize))
    class(ans) <- c("repSparseIdent", class(ans))
    return(ans)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseIdent
setOldClass("repSparseIdent")
setIs("repSparseIdent", "repSparse")

##' Triangular repeated sparse matrix
##'
##' @param diagVals values for the diagonal
##' @param offDiagVals values for the off-diagonal
##' @param low lower triangular?
##' @rdname specialRepSparse
##' @export
repSparseTri <- function(diagVals, offDiagVals, low = TRUE) {
    matSize <- length(diagVals)
    diagIndices <- 1:matSize
    rowIndices <- rep(diagIndices, diagIndices)
    colIndices <- sequence(diagIndices)
    diagIndices <- rowIndices == colIndices
    vals <- numeric(length(diagIndices))
    vals[ diagIndices] <- diagVals
    vals[!diagIndices] <- offDiagVals
    ans <- repSparse(rowIndices, colIndices, seq_along(vals), vals)
    if(!low) ans <- t(ans)
    class(ans) <- c("repSparseTri", class(ans))
    return(ans)
}

##' @name repSparse-class
##' @rdname repSparse-class
##' @family repSparseTopics
##' @exportClass repSparseTri
setOldClass("repSparseTri")
setIs("repSparseTri", "repSparse")

##' Random repeated sparse matrix
##'
##' @param nrows,ncols numbers of rows and columns
##' @param nvals number of values
##' @param nnonzeros number of nonzero elements
##' @param rfunc random number function
##' @param ... dots
##' @rdname specialRepSparse
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
##' @family repSparseTopics
##' @export
chol.repSparse <- function(x, ...) {
    as.repSparse(chol(as.matrix(x, sparse = TRUE)))
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
                     trans = mkCholOneOffDiagTrans(x$vals),
                     Dim = dim(x))
    return(ans)
}

##' @rdname chol
##' @export
##' @examples
##' x <- repSparseCompSymm(1.2, -0.11, 4)
##' as.matrix(chol(x))
##' t(chol(as.matrix(x)))
chol.repSparseCompSymm <- function(x, ...) {
    n <- nrow(x)
    md0 <- x$vals[1]
    od0 <- x$vals[2]
    md <- numeric(n) # main diagonal
    od <- numeric(n-1) # off diagonal
    md[1] <- md0
    od[1] <- od0/sqrt(md0)
    for(i in 2:n) {
        md[i] <- md[i-1] - od[i-1]^2
        if(i == n) break
        od[i] <- (od0 - md0 + md[i])/sqrt(md[i])
    }
    sa <- seq_along(md)
    ii <- rep.int(1:(n-1), (n-1):1)
    ri <- c(sa, unlist(lapply(2:n, ":", n)))
    ci <- c(sa, ii)
    vi <- c(sa, ii + n)
    va <- c(sqrt(md), od)
    repSparse(ri, ci, vi, va, Dim = dim(x))
}

