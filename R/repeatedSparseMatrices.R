## ----------------------------------------------------------------------
## Repeated sparse matrices
##
## used by: strucParseFormula.R, reTrmStructs.R, outputModule.R
## ----------------------------------------------------------------------


## TOC:  grep -A 1 "## \-\-\-" * | grep "\-##" | grep repSparse

##' Repeated sparse matrices
##'
##' A sparse matrix class for matrices whose elements come from a set
##' of a relatively small number of values or computable from a
##' relatively small parameter vector.
##'
##' These parameterized matrix objects are useful for updating
##' matrices at each iteration of a general nonlinear optimizer.  This
##' \code{\link{repSparse-class}} is designed to work well with the
##' \code{\link{Matrix}} package, in that they can be quickly coerced
##' to \code{\link{sparseMatrix}} objects.  The
##' \code{\link{mkSparseTrans}} function can be used to automatically
##' generate a function for updating the values of the
##' \code{sparseMatrix} object from a relatively small number of
##' parameters or repeated values.
##'
##' Such matrices are best constructed by \code{\link{bind}}ing,
##' \code{\link{kron}}ing, and \code{\link{kr}}ing together a series
##' of simpler \code{repSparse} matrices, which can be constructed
##' with the \code{repSparse} function or a function for constructing
##' special \code{repSparse} objects
##' (e.g. \code{\link{repSparseDiag}}).
##'
##' Note that the ordering of the indices supplied to \code{rowInds},
##' \code{colInds}, and \code{valInds} may not appear in the same
##' order in the output.  This is because \code{sortFun} is used to
##' process the order, which by default is \code{\link{standardSort}}.
##' To suppress this behaviour a different \code{sortFun} may be
##' supplied, but this is generally not recommended because certain
##' functions assume certain orderings (e.g. \code{\link{kr}}).
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
##' @param sortFun a function with which to sort the indices of the
##' resulting matrix.  please note that by default the
##' \code{standardSort} function is used, which provides a convenient
##' ordering for computing Khatri-Rao products (with
##' \code{\link{kr}}).  see details.
##' @return A member of the \code{\link{repSparse-class}}.
##' @rdname repSparse
##' @export
##' @examples
##' repSparse(1:4, 4:1, rep(1:2, 2), c(-pi, pi))
repSparse <- function(rowInds, colInds, valInds, vals, trans, Dim,
                      sortFun = standardSort) {
    maxValInds <- ifelse(length(valInds), max(valInds), 0L)
    if(missing(Dim)) Dim <- c(max(rowInds), max(colInds))
    if(!all(seq_len(maxValInds) %in% valInds)) stop("max(valInds) unnecessarily large")
    if(length(vals)    != maxValInds)          stop("mismatch between vals and valInds")
    if(length(rowInds) != length(colInds))     stop("row and column index mismatch")
    if(length(rowInds) != length(valInds))     stop("row and value index mismatch")
    if(missing(trans)) trans <- mkIdentityTrans(vals)
    sortFun(structure(list(rowInds = as.integer(rowInds - 1),
                           colInds = as.integer(colInds - 1),
                           valInds = as.integer(valInds),
                           vals = vals,
                           trans = trans),
                      class = "repSparse",
                      Dim = Dim))
}

setOldClass("repSparse")



## ----------------------------------------------------------------------
## repSparse-class
## ----------------------------------------------------------------------

##' A repeated sparse matrix class
##'
##' An \code{S3} class for repeated sparse matrices, which are lists
##' with a \code{Dim} attribute and with the following elements:
##' \describe{
##'   \item{\code{rowInds}}{0-based row indices}
##'   \item{\code{colInds}}{0-based column indices}
##'   \item{\code{valInds}}{1-based indices for associating each \code{rowInds} and
##'                         \code{colInds} pair with the elements in \code{vals}}
##'   \item{\code{trans}}{A function for transforming a parameter vector
##'                       into \code{vals}.  This function is used by
##'                       \code{\link{update.repSparse}}}
##' }
##' Objects in this class can be constructed with
##' \code{\link{repSparse}} (and, for example,
##' \code{\link{repSparseDiag}}), and has several methods including
##' the following.
##'
##' The \code{...} arguments for the \code{update.repSparse} method
##' can be used to specify special parameter arguments, if
##' \code{object$mkNewPars} exits (which is often the case for special
##' matrices such as \code{\link{repSparseTri}}.  This is a more
##' explicit approach and can therefore be more convenient than having
##' to figure out what order the parameters should appear in
##' \code{newPars}.  For example, \code{update(., diagVals = c(1, 1),
##' offDiagVals = -0.2)} is more explicit than \code{update(., c(1,
##' -0.2, 1))}.
##'
##' @name repSparse-class
##' @rdname repSparse-class
##' @exportClass repSparse
##' @examples
##' (X <- repSparse(rowInds = 1:6,
##'                 colInds = rep(1:2, 3),
##'                 valInds = rep(1:2, each = 3),
##'                 vals    = c(-pi, pi)))
##' (Y <- update(X, c(0.2, 0.8)))
##' as.matrix(X, sparse = TRUE)
##' as.matrix(Y, sparse = TRUE)
##' image(kron(t(X), Y))
##' image(bind(X, Y, type = "diag"))
setOldClass("repSparse")

##' @param x \code{repSparse} object
##' @param n how many indices and repeated values to print?
##' @param ... passed to subsequent functions
##' @rdname repSparse-class
##' @method print repSparse
##' @export
print.repSparse <- function(x, n = 6L, ...) {

    init <- getInit(x)
    parsEqualRepVals <- identical(x$vals, init)
    headVals <- head(x$vals, n = n)
    inds <- as.data.frame(x)
    headInds <- head(inds, n = n)
    reportTruncVals <- length(x$vals) > n
    reportTruncInds <- nrow(inds) > n
    reportTruncPars <- length(init) > n

    cat("Repeated sparse matrix\n")
    cat("----------------------\n")
    cat("dimensions:\n")
    cat("", nrow(x), "rows and", ncol(x), "columns\n",
        nrow(inds), "nonzero values\n",
        length(x$vals), "repeated values\n")
    if(!parsEqualRepVals) {
        cat("", length(init), "parameters\n")
    }
    cat("\n")
    if(reportTruncVals) cat("first", n, "")
    cat("repeated values:\n")
    print(headVals)
    if(!parsEqualRepVals) {
        cat("\n")
        if(reportTruncPars) cat("first", n, "")
        cat("initial parameters:\n")
        print(head(getInit(x), n))
    }
    cat("\n")
    if(reportTruncInds) cat("first", n, "")
    cat ("row, column, and value indices (with associated values):\n")
    print(headInds)
    if(length(class(x)) > 1L) {
        cat("\nspecial repeated sparse matrix, inheriting from:\n")
        print(class(x))
    }
    cat("\n")
}

##' @param object \code{repSparse} object
##' @param newPars new parameter values
##' @importFrom stats update
##' @rdname repSparse-class
##' @method update repSparse
##' @export
update.repSparse <- function(object, newPars, ...) {
    l... <- list(...)
    if(length(l...) > 0L) {
        mkNewPars <- object$mkNewPars
        if(is.null(mkNewPars)) stop("special arguments given, ",
                                    "but no function available for ",
                                    "constructing a parameter vector")
        newPars <- do.call(mkNewPars, l...)
    }
    if(missing(newPars)) {
        if(length(gi <- getInit(object)) == 0L) {
            return(object)
        }
        newPars <- gi
    }
    if(length(newPars) != (np <- length(getInit(object))))
        stop("newPars must have the same length as getInit(object), ",
             "which in this case is ", np)
    object$vals <- object$trans(newPars)
    return(object)
}

.removeLatticeWhitespace <- function() {
    lh <- lattice::trellis.par.get("layout.heights")
    lw <- lattice::trellis.par.get("layout.widths")
    lh[grep("padding", names(lh))] <- lw[grep("padding", names(lw))] <- 0
    ac <- lattice::trellis.par.get("axis.components")
    for(i in 1:4) ac[[i]][c("pad1", "pad2")] <- 0
    lattice::trellis.par.set(layout.heights = lh, layout.widths = lw,
                              axis.components = ac)
}

.xscaleComponents <- function(...) {
    ans <- lattice::xscale.components.default(...)
    ans$bottom$labels$labels <- rep("", length(ans$bottom$labels$labels))
    ans$bottom$ticks$tck <- 0
    #ans$bottom$ticks$at <- c()
    #ans$bottom$labels$check.overlap <- FALSE
    #ans$top <- FALSE
    ans
}

.yscaleComponents <- function(...) {
    ans <- lattice::yscale.components.default(...)
    ans$left$labels$labels <- rep("", length(ans$left$labels$labels))
    ans$left$ticks$tck <- 0
    #ans$left$ticks$at <- c()
    #ans$left$labels$check.overlap <- FALSE
    #ans$right <- FALSE
    ans
}

##' @param plain should a completely plain plot be used? (try and see)
##' @importFrom Matrix image
##' @rdname repSparse-class
##' @method image repSparse
##' @export
image.repSparse <- function(x, plain = FALSE, ...) {
    blank <- length(x$vals) == 0
    x <- as.matrix(x, sparse = TRUE)
    if(plain) {
        if(blank) { # hack to give something useful for
                                  # image(repSparseBlank(...))
            grid::pushViewport(grid::viewport())
            bx <- grid::unit(0.99, "npc")
            grid::grid.rect(width = bx, height = bx)
        } else {
            .removeLatticeWhitespace()
            image(x, sub = "", xlab = "", ylab = "", colorkey = FALSE,
                  xscale.components = .xscaleComponents,
                  yscale.components = .yscaleComponents)
        }
    } else {
        if(blank) {
            stop("can't plot blank matrix in this way, ",
                 "try image(..., plain = TRUE)")
        } else {
        ## not sure why i can't pass ... directly, but apparently this
        ## doesn't work
            do.call(image, c(list(x), list(...)))
        }
    }
}

##' @param y not used
##' @rdname repSparse-class
##' @method plot repSparse
##' @export
plot.repSparse <- function(x, y, plain = FALSE, ...) image(x, plain = plain, ...)

##' @rdname repSparse-class
##' @method t repSparse
##' @export
t.repSparse <- function(x) {
    tx <- x
    tx$rowInds <- x$colInds
    tx$colInds <- x$rowInds
    attr(tx, "Dim") <- rev(dim(x))
    return(standardSort(tx))
}



##' @rdname repSparse-class
##' @export
setMethod("diag", signature(x = "repSparse"), {
    function(x) {
        with(x, {
            diagInds <- which(rowInds == colInds)
            rowDiagInds <- rowInds[diagInds] + 1L
            fullInds <- match(seq_len(min(dim(x))), rowDiagInds)
            ans <- vals[valInds[diagInds]][fullInds]
            ifelse(is.na(ans), 0, ans)
        })
    }
})

##' @param value replacement value for diagonal
##' @rdname repSparse-class
##' @export
setMethod("diag<-", signature(x = "repSparse"), {
    function(x, value) {
        diagInds <- with(x, {
            diagInds <- unique(valInds[  rowInds == colInds])
            offDiagInds <- unique(valInds[!(rowInds == colInds)])
            if(length(diagInds) != length(value)) { # length mismatch
                stop("too many replacement values for the number of repeated diagonal values")
            }
            if(any(diagInds %in% offDiagInds)) { # complex pattern
                stop(" the repetition pattern is too complicated to replace the diagonal.\n",
                     "some repeated values are both on the diagonal and off of it")
            }
            diagInds
        })
        x$vals[diagInds] <- value
        return(x)
    }
})


##' @rdname repSparse-class
##' @method dim repSparse
##' @export
dim.repSparse <- function(x) attr(x, "Dim")


##' Sort the indices of a repeated sparse matrix
##'
##' The ordering of the indices is arbitrary, but some orders are more
##' computationally efficient than others in certain circumstances.
##' \code{standardSort} is defined as \code{sort(sort(., type =
##' "row"), type = "col")}, which is often the 'best' order.  It is
##' used throughout, and in particular is required for using the
##' \code{\link{kr}} function.
##'
##' @param x \code{\link{repSparse}} object
##' @param decreasing see \code{\link{sort}}
##' @param type sort by column, row, or value indices?
##' @param ... for consistency with generic
##' @rdname sort.repSparse
##' @method sort repSparse
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

##' @rdname sort.repSparse
##' @export
standardSort <- function(x) {
    sort(sort(x, type = "row"), type = "col")
}


## ----------------------------------------------------------------------
## Coercion -- as...
## ----------------------------------------------------------------------

##' Coerce to and from repeated sparse matrices
##'
##' @param x an object
##' @param ... dots
##' @rdname as.repSparse
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
as.repSparse.dist <- function(x, ...) {
    matSize <- attr(x, "Size")
    lowTriSize <- matSize - 1
    diagVals <- rep(0, matSize)
    offDiagInds <- triInds(rep(1:lowTriSize, 1:lowTriSize),
                           sequence(1:lowTriSize),
                           lowTriSize)
    repSparseSymm(diagVals, x[offDiagInds])
}

##' @rdname as.repSparse
##' @export
as.repSparse.dCHMsimpl <- function(x, ...) {
    matSize <- x@Dim[1]
    rowInds <- rep(1:matSize, 1:matSize)
    colInds <- sequence(1:matSize)
    valInds <- triInds(rowInds, colInds, matSize)
    vals <- as(x, "sparseMatrix")@x
    repSparse(rowInds, colInds, valInds, vals)
}


##' @rdname as.repSparse
##' @method as.data.frame repSparse
##' @export
as.data.frame.repSparse <- function(x, ...) {
    with(x, {
        data.frame(rowInds = rowInds,
                   colInds = colInds,
                   valInds = valInds,
                   vals = vals[valInds])
    })
}


##' @importFrom Matrix sparseMatrix
##' @rdname as.repSparse
##' @param sparse return \code{sparseMatrix}?
##' @method as.matrix repSparse
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

##' @rdname as.repSparse
##' @method as.matrix repSparseChol
##' @export
as.matrix.repSparseChol <- function(x, sparse = FALSE, ...) {
    ans <- as.matrix.repSparse(x, sparse = sparse, ...)
    return(as(ans, "dtCMatrix"))
}


repSparse2gCsparse <- function(from) {
    ans <- new("dgCMatrix")
    ans@i <- as.integer(from$rowInds)
    ans@p <- as.integer(ind2point(from$colInds, ncol(from)))
    ans@x <- as.numeric(with(from, vals[valInds]))
    ans@Dim <- as.integer(dim(from))
    return(ans)
}

repSparse2gTsparse <- function(from) {
    ans <- new("dgTMatrix")
    ans@i <- as.integer(from$rowInds)
    ans@j <- as.integer(from$colInds)
    ans@x <- as.numeric(with(from, vals[valInds]))
    ans@Dim <- as.integer(dim(from))
    return(ans)
}

repSparse2tCsparse <- function(from) {
    as(as(from, "dgCMatrix"), "dtCMatrix")
}

repSparse2tTsparse <- function(from) {
    as(as(from, "dgTMatrix"), "dtTMatrix")
}

##' as("repSparse", "sparseMatrix")
##' @name repSparse2sparseMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix sparseMatrix
##' @importFrom methods coerce
##' @importFrom methods as
setAs("repSparse",  "sparseMatrix", def = repSparse2gCsparse)

##' as("repSparse", "CsparseMatrix")
##' @name repSparse2gCsparseMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix CsparseMatrix
setAs("repSparse", "CsparseMatrix", def = repSparse2gCsparse)

##' as("repSparse", "dgCMatrix")
##' @name repSparse2dgCMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix dgCMatrix
setAs("repSparse",     "dgCMatrix", def = repSparse2gCsparse)

##' as("repSparse", "dtCMatrix")
##' @name repSparse2dtCMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix dtCMatrix
setAs("repSparse",     "dtCMatrix", def = repSparse2tCsparse)

##' as("repSparse", "TsparseMatrix")
##' @name repSparse2gTsparseMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix TsparseMatrix
setAs("repSparse", "TsparseMatrix", def = repSparse2gTsparse)

##' as("repSparse", "dgTMatrix")
##' @name repSparse2dgTMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix dgTMatrix
setAs("repSparse", "dgTMatrix", def = repSparse2gTsparse)

##' as("repSparse", "dtTMatrix")
##' @name repSparse2dtTMatrix
##' @rdname as.repSparse
##' @importClassesFrom Matrix dtTMatrix
setAs("repSparse",     "dtTMatrix", def = repSparse2tTsparse)



##' Get repetition pattern of a repeated sparse matrix
##'
##' @param object a \code{\link{repSparse}} object
##' @export
getRepPattern <- function(object) {
    object$vals <- seq_along(object$vals)
    return(as.matrix(object, sparse = TRUE))
}




## ----------------------------------------------------------------------
## Matrix operations -- kron (Kronecker product), kr (Khatri-Rao
## product)
## ----------------------------------------------------------------------

##' Kronecker and Khatri-Rao products for repeated sparse matrices
##'
##' @param X,Y repeated sparse matrices (\code{\link{repSparse-class}})
##' @param trans see argument \code{FUN} in \code{\link{outer}}
##' @param makedimnames ignored
##' @param ... ignored
##' @rdname matrixOperations
##' @family matrixCombining
##' @export
##' @examples
##' set.seed(1)
##' X <- repSparse(c(1, 2, 1, 2), c(1, 1, 2, 2), 1:4, rnorm(4))
##' Y <- repSparse(c(1, 2, 1), c(1, 1, 2), 1:3, rnorm(3))
##' (kronExample <- kron(X, Y))
kron <- function(X, Y, trans = "*",
                 makedimnames = FALSE,  ...) {

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
              Dim = dim(X) * dim(Y))
}

setOldClass("repSparseKron")
setIs("repSparseKron", "repSparse")


##' @rdname matrixOperations
##' @family matrixCombining
##' @export
##' @examples
##' (krExample <- kr(X, Y))
kr <- function(X, Y, trans = "*") {

    ## modified Matrix::KhatriRao to allow for repeated sparse case

    X <- standardSort(X)
    Y <- standardSort(Y)

    p <- ncol(X)
    xp <- as.integer(ind2point(X$colInds, p))
    yp <- as.integer(ind2point(Y$colInds, ncol(Y)))
    xn <- diff(xp)
    yn <- diff(yp)
    newp <- as.integer(diffinv(xn * yn))

    xn.yp <- xn[as.logical(yn)]
    yj <- factor(Y$colInds)
    rep.yn <- rep.int(yn, xn)
    i1 <- rep.int(X$rowInds, rep.yn)
    i2 <- unlist(rep(split.default(Y$rowInds, yj), xn.yp))
    n1 <- nrow(X); n2 <- nrow(Y)
    dim <- as.integer(c(n1 * n2, p))

    v1 <- rep.int(X$valInds, rep.yn)
    v2 <- unlist(rep(split.default(Y$valInds, yj), xn.yp))

    newRowInds <- i1 * n2 + i2
    newColInds <- point2ind(newp) - 1L
    newValInds <- (v1 - 1L) * length(Y$vals) + v2
    newVals <- as.vector(outer(Y$vals, X$vals, FUN = trans))
    newTrans <- mkOuterTrans(Y$trans, X$trans, trans)

    structure(list(rowInds = newRowInds,
                   colInds = newColInds,
                   valInds = newValInds,
                   vals = newVals,
                   trans = newTrans),
              class = c("repSparseKr", "repSparse"),
              Dim = c(dim(X)[1] * dim(Y)[1], dim(X)[2]))
}

setOldClass("repSparseKr")
setIs("repSparseKr", "repSparse")

##' @rdname matrixOperations
##' @export
setMethod("kronecker", signature(X = "repSparse", Y = "repSparse"), {
    function(X, Y, FUN = "*", make.dimnames = FALSE, ...) kron(X, Y)
})

##' @param e1,e2 \code{repSparse} objects
##' @note The \code{*} operator is matrix multiplication, not
##' element-wise multiplication.
##' @rdname matrixOperations
##' @export
Ops.repSparse <- function(e1, e2) {
    FUN = .Generic
    if(FUN == "*") {
        return(as.matrix(e1, sparse = TRUE) %*% as.matrix(e2, sparse = TRUE))
    } else {
        return(getGeneric(FUN)(as.matrix(e1, sparse = TRUE), as.matrix(e2, sparse = TRUE)))
    }
}

##' @param x,y \code{\link{repSparse}} matrix objects
##' @rdname matrixOperations
##' @importFrom Matrix crossprod
##' @export
setMethod("crossprod", signature(x = "repSparse", y = "repSparse"), {
    function(x, y = NULL, ...) {
        crossprod(as.matrix(x, sparse = TRUE),
                  as.matrix(y, sparse = TRUE), ...)
    }
})

##' @param x,y \code{\link{repSparse}} matrix objects
##' @rdname matrixOperations
##' @importFrom Matrix tcrossprod
##' @export
setMethod("tcrossprod", signature(x = "repSparse", y = "repSparse"), {
    function(x, y = NULL, ...) {
        tcrossprod(as.matrix(x, sparse = TRUE),
                   as.matrix(y, sparse = TRUE), ...)
    }
})

## ----------------------------------------------------------------------
## Make trans functions
## ----------------------------------------------------------------------

##' Construct functions for transforming a parameter vector to the
##' non-zero values of a repeated sparse matrix
##'
##' These functions return a 'trans function', for transforming a
##' parameter vector to the repeated non-zero values of a repeated
##' sparse matrix.  Each trans function takes one vector argument
##' called \code{matPars}, which are the parameters of the matrix.
##' The environment of these transformation functions must contain a
##' vector called \code{init}, which contains the initial values (aka,
##' the prototype) of \code{matPars}.
##'
##' @rdname mkTrans
##' @aliases mkTrans
##' @export
mkIdentityTrans <- function(init) {
    local({
        init <- init
        function(matPars) matPars
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
        n <- nChoose2Inv(m) - 1 # as.integer((sqrt(1 + 8 * m) - 1)/2)
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


##' @rdname mkTrans
##' @export
mkConstVarCholTrans <- function(init) {
    local({
        init <- init
        matSize <- nChoose2Inv(length(init) - 1L) # minus one for the
                                                  # standard deviation
                                                  # parameter
        diagIndices <- 1:matSize
        rowIndices <- rep(diagIndices, diagIndices)
        colIndices <- sequence(diagIndices)
        diagIndices <- rowIndices == colIndices
        vals <- numeric(length(diagIndices))
        function(matPars) {
            sdVal <- matPars[1]
            offDiagVals <- matPars[-1]
            innProd <- (sdVal^2) -
                c(0, tapply(offDiagVals^2, rowIndices[!diagIndices], sum))
            if(any(innProd < 0L)) {
                stop("resulting matrix not positive definite")
            }
            diagVals <- sqrt(innProd)
            vals[ diagIndices] <- diagVals
            vals[!diagIndices] <- offDiagVals
            return(vals)
        }
    })
}

##' @rdname mkTrans
##' @export
mkCorMatCholTrans <- function(init) {
    local({
        init <- init
        matSize <- nChoose2Inv(length(init))
        diagIndices <- 1:matSize
        rowIndices <- rep(diagIndices, diagIndices)
        colIndices <- sequence(diagIndices)
        diagIndices <- rowIndices == colIndices
        vals <- numeric(length(diagIndices))
        sdVal <- 1
        offDiagFun <- function(x) {
            x <- x^2
            sqrt(x / (sum(x) + 1))
        }
        innProdFun <- function(x) sum(x^2)
        function(matPars) {
            splitPars <- split(matPars, rowIndices[!diagIndices])
            offDiagValsList <- lapply(splitPars, offDiagFun)
            innProd <- 1 - c(0, sapply(offDiagValsList, innProdFun))
            if(any(innProd < 0L)) stop("resulting matrix not positive definite") ## shouldn't ever happen
            diagVals <- sqrt(innProd)
            ## a little DRY
            vals[ diagIndices] <- diagVals
            vals[!diagIndices] <- sign(matPars) * unlist(offDiagValsList)
            return(vals)
        }
    })
}

##' @param cholObj,symmObj corresponding triangular cholesky and
##' symmetric \code{\link{repSparse}} matrix objects
##' @param vecDist distances as a vector
##' @rdname mkTrans
##' @export
mkExpCholTrans <- function(init, cholObj, symmObj, vecDist) {
    local({
        init <- init
        cholObj <- cholObj
        symmObj <- symmObj
        vecDist <- as.numeric(vecDist) + 0
        symmObj@x <- exp(-init * vecDist)
        function(matPars) {
            symmObj@x <- exp(-matPars * vecDist)
            as(update(cholObj, symmObj), "sparseMatrix")@x
        }
    })
}


##' @rdname mkTrans
##' @export
##' @param defaultOutput default vector for diagonal elements
##' @param devfunEnv environment of the deviance function
mkFlexDiagTrans <- function(init, defaultOutput,
                            devfunEnv) {
    local({
        init <- init
        nBasis <- length(init)
        if(nBasis == 0L) stop("initial parameter vector must have length greater than zero")
        defaultOutput <- defaultOutput
        devfunEnv <- devfunEnv
        Xspline <- NULL
        if(nBasis == 1) {
            function(matPars) {
                return(as.numeric(rep(exp(matPars), length(defaultOutput))))
            }
        } else {
            function(matPars) {
                lp <- try(evalq(pp$linPred(1) + resp$offset, devfunEnv), silent = TRUE)
                b1 <- try(evalq(pp$b      (1),               devfunEnv), silent = TRUE)
                if(inherits(lp, "try-error")) return(defaultOutput)
                lp <- try(lp - b1[indsObsLevel], silent = TRUE)
                if(inherits(lp, "try-error")) return(defaultOutput)
                Xspline <- try(ns(lp, nBasis, intercept = TRUE), silent = TRUE)
                assign("Xspline", Xspline, envir = parent.env(environment()))
                if(inherits(Xspline, "try-error")) return(defaultOutput)
                return(as.numeric(exp(Xspline %*% matPars)))
            }
        }
    })
}


##' @param Atrans,Btrans functions for transforming two repeated
##' sparse matrices, \code{A} and \code{B} say
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

##' @param transList list of transformation functions
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


##' @param covariate covariate
##' @param grpFac a grouping factor (or anything coercible to
##' \code{numeric} really) with the number of levels equal to the
##' length of \code{init}
##' @rdname mkTrans
##' @export
mkVarExpTrans <- function(init, covariate, grpFac) {
    local({
        init <- init
        covariate <- covariate
        grpFac <- grpFac
        function(matPars) {
            parsExpand <- matPars[as.numeric(grpFac)]
            return(exp(2 * parsExpand * covariate))
        }
    })
}

##' @param modMat model matrix
##' @rdname mkTrans
##' @export
mkMultiVarExpTrans <- function(init, modMat, grpFac) {
    local({
        init <- init
        if(!is.na(grpFac)) {
            termMat <- t(KhatriRao(as(grpFac, "sparseMatrix"),
                                   t(as(modMat, "sparseMatrix"))))
        } else {
            termMat <- modMat
        }
        function(matPars) {
            return(exp(2 * as.numeric(termMat %*% matPars)))
        }
    })
}

##' @rdname mkTrans
##' @export
mkVarPowerTrans <- function(init, covariate, grpFac) {
    local({
        init <- init
        covariate <- covariate
        grpFac <- grpFac
        function(matPars) {
            parsExpand <- matPars[as.numeric(grpFac)]
            return(abs(covariate)^(2 * parsExpand))
        }
    })
}


##' @rdname mkTrans
##' @export
mkVarIdentTrans <- function(init, covariate, grpFac) {
    local({
        init <- init
        grpFac <- grpFac
        function(matPars) {
            return(matPars[as.numeric(grpFac)])
        }
    })
}


##' @rdname mkTrans
##' @export
mkVarFixedTrans <- function(init, covariate, grpFac) {
    local({
        init <- init
        covariate <- covariate
        function(matPars) {
            return(abs(covariate))
        }
    })
}


##' @param matSize number of rows (or columns) of the matrix
##' @rdname mkTrans
##' @export
mkCholCompSymmTrans <- function(init, matSize) {
    local({
        init <- init
        matSize <- matSize
        function(matPars) {
            md0 <- matPars[1]
            od0 <- matPars[2]
            md <- numeric(matSize) # main diagonal
            od <- numeric(matSize-1) # off diagonal
            md[1] <- md0
            od[1] <- od0/sqrt(md0)
            for(i in 2:matSize) {
                md[i] <- md[i-1] - od[i-1]^2
                if(i == matSize) break
                od[i] <- (od0 - md0 + md[i])/sqrt(md[i])
            }
            return(c(sqrt(md), od))
        }
    })
}


## ----------------------------------------------------------------------
## Special modifications of repeated sparse matrices
## ----------------------------------------------------------------------

##' Constant repeated sparse matrix
##'
##' Convert a parameterized repeated sparse matrix into a constant
##' repeated sparse matrix (i.e. \code{getInit(object)} has length
##' zero).
##'
##' @param object repeated sparse matrix
##' @family modifications
##' @export
resetTransConst <- function(object) {
    object$trans <- local({
        baseline <- object$vals
        init <- numeric(0)
        function(matPars = NULL) {
            return(baseline)
        }
    })
    ## FIXME: return(simplifyRepSparse(object)) ???
    return(object)
}


##' Simplify repeated sparse matrices
##'
##' Simplify repeated sparse matrices by (1) forcing non structural
##' zeros to be structural zeros, (2) eliminate duplicates in the
##' repeated values, and (3) resetting the transformation function to
##' be the identity.
##'
##' @param object repeated sparse matrix
##' @param force force non structural zeros to be structural zeros?
##' @param eliminate eliminate duplicates in the repeated values?
##' @param reset reset the transformation function to be the identity?
##' @param ... not yet used
##' @family modifications
##' @export
simplifyRepSparse <- function(object,
                              force = TRUE,
                              eliminate = TRUE,
                              reset = TRUE, ...) {
                                        # force non structural zeros
                                        # to be structural
    if(force) {
        valsToKeep <- which(object$vals != 0L)
        elementsToKeep <- with(object, valInds %in% valsToKeep)
        object$vals <- object$vals[valsToKeep]
        object$rowInds <- object$rowInds[elementsToKeep]
        object$colInds <- object$colInds[elementsToKeep]
        object$valInds <- flattenIntVec(object$valInds[elementsToKeep])
    }
                                        # eliminate duplicates in the
                                        # repeated values
    if(eliminate) {
        newVals <- with(object, vals[which(!duplicated(vals))])
        newValInds <- with(object, match(vals[valInds], newVals))
        object$vals <- newVals
        object$valInds <- newValInds
    }
                                        # reset to identity trans
    if(reset) {
        object$trans <- mkIdentityTrans(object$vals)
    }
    return(object)
}

##' Add scalar multiple to a repeated sparse matrix
##'
##' Adjust the \code{trans} function of \code{object} such that a
##' scalar multiplier parameter is concatenated at the beginning of
##' the parameter vector.
##'
##' @param object repeated sparse matrix
##' @param mult initial value for the multiplier parameter
##' @family modifications
##' @export
scalarMult <- function(object, mult) {
    object$trans <- local({
        init <- c(mult, getInit(object))
        unscaledTrans <- object$trans
        function(matPars) {
            matPars[1] * unscaledTrans(matPars[-1])
        }
    })
    ans <- update(object)
    class(ans) <- c("repSparseScalarMult", class(ans))
    return(ans)
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
##' @rdname bind
##' @family matrixCombining
##' @export
##' @examples
##' set.seed(1)
##' X <- repSparse(c(1, 2, 1, 2), c(1, 1, 2, 2), 1:4, rnorm(4))
##' Y <- repSparse(c(1, 2, 1), c(1, 1, 2), 1:3, rnorm(3))
##' Z <- repSparse(c(1, 2), c(1, 2), 1:2, rnorm(2))
##' (bindDiag <- bind(X, Y, Z, type = "diag"))
##' (bindRow  <- bind(X, Y, Z, type = "row"))
##' (bindCol  <- bind(X, Y, Z, type = "col"))
bind <- function(..., type = c("row", "col", "diag")) {
    mats <- list(...)
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
        return(ans)
    })
}

setOldClass("repSparseBind")
setIs("repSparseBind", "repSparse")


##' @param mats list of \code{repSparse} matrix objects
##' @rdname bind
##' @export
.bind <- function(mats, type = c("row", "col", "diag")) {
    do.call(bind, c(mats, list(type = type)))
}

##' Repeat repeated sparse matrix
##'
##' @param x \code{repSparse} object
##' @param times like \code{rep}
##' @param repVals should the unique values (\code{x$vals}) be
##' repeated too?  in other words, do you want to have each replicate
##' to have guaranteed identical values (\code{FALSE}, the default),
##' or do you want each replicate to be settable to it's own values
##' (\code{TRUE})?
##' @rdname bind
##' @family matrixCombining
##' @export
##' @examples
##' set.seed(1)
##' X <- repSparse(c(1, 2, 1, 2), c(1, 1, 2, 2), 1:4, rnorm(4))
##' Y <- repSparse(c(1, 2, 1), c(1, 1, 2), 1:3, rnorm(3))
##' Z <- repSparse(c(1, 2), c(1, 2), 1:2, rnorm(2))
##' (repDiag <- rep(X, 3, type = "diag"))
##' (repRow  <- rep(X, 3, type = "row"))
##' (repCol  <- rep(X, 3, type = "col"))
rep.repSparse <- function(x, times,
                          type = c("row", "col", "diag"),
                          repVals = FALSE,
                          ...) {

    type = type[[1]]
    len <- length(x$rowInds)
    rowInds <- rep.int(x$rowInds, times)
    colInds <- rep.int(x$colInds, times)
    rowMult <- colMult <- 1
    if((type == "row") || (type == "diag")) {
        off <- rep.int((0:(times - 1)) * nrow(x), rep.int(len, times))
        rowInds <- rowInds + off
        rowMult <- times
    }
    if((type == "col") || (type == "diag")) {
        off <- rep.int((0:(times - 1)) * ncol(x), rep.int(len, times))
        colInds <- colInds + off
        colMult <- times
    }
    if(repVals) {
        lenXvals <- length(x$vals)
        vals <- rep(x$vals, times)
        off <- rep.int((0:(times - 1)) * lenXvals,
                       rep.int(len, times))
        valInds <- x$valInds + off
        trans <- mkListTrans(rep(list(x$trans), times))
    } else {
        vals <- x$vals
        valInds <- rep(x$valInds, times = times)
        trans <- x$trans
    }
    ans <- repSparse(rowInds = rowInds + 1,
                     colInds = colInds + 1,
                     valInds = valInds,
                     vals = vals,
                     trans = trans,
                     Dim = c(rowMult, colMult) * dim(x))
    class(ans) <- c("repSparseRep", class(ans))
    ans$mkNewPars <- x$mkNewPars
    return(ans)
}

setOldClass("repSparseRep")
setIs("repSparseRep", "repSparse")


## ----------------------------------------------------------------------
## Subsetting
## ----------------------------------------------------------------------

##' Subsetting repeated sparse matices
##'
##' @param x a \code{\link{repSparse-class}} object
##' @param rowInds,colInds 1-based integer indices for rows and
##' columns
##' @param ... unused
##' @family matrixCombining
##' @rdname subset
##' @export
##' @examples
##' set.seed(1)
##' X <- repSparse(c(1, 2, 1, 2), c(1, 1, 2, 2), 1:4, rnorm(4))
##' (subsetExample <- subset(X, c(1, 2, 2, 2, 1, 1)))
subset.repSparse <- function(x, rowInds = NULL, colInds = NULL, ...) {
    if(!is.null(rowInds)) x <-   repSparseRowSubset(  x , rowInds)
    if(!is.null(colInds)) x <- t(repSparseRowSubset(t(x), colInds))
    return(x)
}


##' @rdname subset
##' @export
repSparseRowSubset <- function(x, rowInds) {

                                        # save the original indices
    ri <- x$rowInds + 1L
    ci <- x$colInds + 1L
    vi <- x$valInds

                                        # find out about the indices
    coir <- countInRange(ri)
    rowIndsCounts <- coir$counts[rowInds]
    tripletsToKeep <- rowInds[rowInds %in% which(coir$counts > 0L)]
    riFac <- factor(ri, levels = coir$vals)

                                        # modify the indices
    x$rowInds <- as.integer(rep(seq_along(rowIndsCounts), rowIndsCounts) - 1L)
    x$colInds <- as.integer(unlist(split(ci, riFac)[tripletsToKeep]) - 1L)
    x$valInds <- as.integer(unlist(split(vi, riFac)[tripletsToKeep]))
    attr(x, "Dim")[1] <- as.integer(length(rowInds))

    return(x)

    ## FIXME: decided not to use this code to remove any possible
    ## unused repeated values, mostly because it hurts my head to
    ## think about how to deal with special parameterized matrics

    ## valsToKeep <- seq_along(va) %in% viNew
    ## vaNew <- va[valsToKeep]
    ## viNew <- flattenIntVec(viNew)

    ## repSparse(riNew, ciNew, viNew, vaNew,
    ##           Dim = c(length(rowInds), ncol(x)),
    ##           trans = x$trans)

}



## ----------------------------------------------------------------------
## Changing sparse formats
## ----------------------------------------------------------------------

##' Changing sparse format
##'
##' \code{point2ind} takes a vector of column pointers, as used in
##' sparse matrices of class \code{dgCMatrix}, and returns a vector of
##' column indices, as used in sparse matrices of class
##' \code{dgTMatrix}.  \code{ind2point} is the approximate inverse of
##' \code{point2ind}.  See \code{\link{sparseMatrix}} for more
##' information about these conversions.
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

##' Blank repeated sparse matrix
##'
##' @param nrow,ncol number of rows and columns
##' @export
##' @family repSparseSpecial
##' @examples
##' (xBlank <- repSparseBlank(5, 5))
repSparseBlank <- function(nrow, ncol) {
    ## MATNAME: Blank
    ans <- repSparse(integer(0), integer(0), integer(0), numeric(0), Dim = c(nrow, ncol))
    class(ans) <- c("repSparseBlank", class(ans))
    return(ans)
}

setOldClass("repSparseBlank")
setIs("repSparseBlank", "repSparse")


##' Identity repeated sparse matrix
##'
##' @param matSize matrix size
##' @family repSparseSpecial
##' @export
##' @examples
##' (xIdent <- repSparseIdent(5))
repSparseIdent <- function(matSize) {
    ## MATNAME: Identity
    ans <- repSparseDiag(1, rep(1, matSize))
    class(ans) <- c("repSparseIdent", class(ans))
    return(ans)
}

setOldClass("repSparseIdent")
setIs("repSparseIdent", "repSparse")


##' Diagonal repeated sparse matrix
##'
##' @param vals vector of values
##' @param valInds vector of value indices
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xDiag <- repSparseDiag(rnorm(5)))
repSparseDiag <- function(vals, valInds) {
    ## MATNAME: Diagonal
    if(missing(valInds)) valInds <- seq_along(vals)
    matSize <- length(valInds)
    ans <- repSparse(1:matSize, 1:matSize, valInds, vals)
    class(ans) <- c("repSparseDiag", class(ans))
    return(ans)
}

setOldClass("repSparseDiag")
setIs("repSparseDiag", "repSparse")



##' Full repeated sparse matrix
##' 
##' @param nrow,ncol numbers of rows and columns
##' @param vals vector of values
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xFull <- repSparseFull(5, 5, rnorm(25)))
repSparseFull <- function(nrow, ncol, vals) {
    ## MATNAME: Full
    ans <- repSparse(rep(1:nrow, ncol),
                     rep(1:ncol, each = nrow),
                     1:(nrow * ncol),
                     vals)
    class(ans) <- c("repSparseFull", class(ans))
    return(ans)       
}

setOldClass("repSparseFull")
setIs("repSparseFull", "repSparse")


##' Column vector as repeated sparse matrix
##'
##' @param vals vector of values
##' @param valInds vector of value indices
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xCol <- repSparseCol(rnorm(5)))
repSparseCol <- function(vals, valInds) {
    ## MATNAME: Column vector
    if(missing(valInds)) valInds <- seq_along(vals)
    matSize <- length(valInds)
    ans <- repSparse(1:matSize, rep(1, matSize), valInds, vals)
    class(ans) <- c("repSparseCol", class(ans))
    return(ans)
}

setOldClass("repSparseCol")
setIs("repSparseCol", "repSparse")


##' Repeated sparse indicator matrix
##'
##' @param fac vector coercible to factor
##' @family repSparseSpecial
##' @export
##' @examples
##' (xInd <- repSparseInd(rep(1:3, c(2, 2, 1))))
repSparseInd <- function(fac) {
    ## MATNAME: Indicator
    ans <- as.repSparse(as.factor(fac))
    class(ans) <- c("repSparseInd", class(ans))
    return(ans)
}

setOldClass("repSparseInd")
setIs("repSparseInd", "repSparse")


##' Triangular repeated sparse matrix
##'
##' @param diagVals values for the diagonal
##' @param offDiagVals values for the off-diagonal
##' @param low lower triangular?
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xTri <- repSparseTri(rnorm(5), rnorm(choose(5, 2))))
repSparseTri <- function(diagVals, offDiagVals, low = TRUE) {
    ## MATNAME: Triangular
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
    ans$mkNewPars <- local({
        diagIndices
        function(diagVals, offDiagVals) {
            vals <- numeric(length(diagIndices))
            vals[ diagIndices] <- diagVals
            vals[!diagIndices] <- offDiagVals
            return(vals)
        }
    })
    return(ans)
}

setOldClass("repSparseTri")
setIs("repSparseTri", "repSparse")

##' General and full triangular repeated sparse matrix
##'
##' @param nrow,ncol
##' @param vals values for the nonzero values
##' @param diag include diagonal?
##' @param low lower triangular?
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xTri <- repSparseGenFullTri(rnorm(5), rnorm(choose(5, 2))))
repSparseGenFullTri <- function(nrow, ncol, vals, diag = TRUE, low = TRUE) {
    ## MATNAME: General full triangular
    rowIndices <- rev(nrow - sequence(nrow - (ncol:1) + 1)) + 1
    colIndices <- rep(0:(ncol - 1), nrow:(nrow - ncol + 1)) + 1
    ans <- repSparse(rowIndices, colIndices, seq_along(vals), vals)
    class(ans) <- c("repSparseGenFullTri", class(ans))
    return(ans)
}

setOldClass("repSparseGenFullTri")
setIs("repSparseGenFullTri", "repSparse")


##' Repeated sparse matrix of ones
##'
##' @param nrow,ncol numbers of rows and columns
##' @family repSparseSpecial
##' @export
##' @examples
##' (xOnes <- repSparseOnes(5, 5))
repSparseOnes <- function(nrow, ncol) {
    ## MATNAME: Ones
    rc <- expand.grid(seq_len(nrow), seq_len(ncol))
    vi <- rep(1, nrow(rc))
    ans <- repSparse(rc[, 1], rc[, 2], vi, 1)
    class(ans) <- c("repSparseOnes", class(ans))
    return(ans)
}

setOldClass("repSparseOnes")
setIs("repSparseOnes", "repSparse")


##' Symmetric repeated sparse matrix
##'
##' @param diagVals values for the diagonal
##' @param offDiagVals values for the off-diagonal
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xSymm <- repSparseSymm(rnorm(5), rnorm(choose(5, 2))))
repSparseSymm <- function(diagVals, offDiagVals) {
    ## MATNAME: Symmetric
    matSize <- length(diagVals)
    diagIndices <- 1:matSize
    repIndices <- rep(diagIndices, diagIndices)
    seqIndices <- sequence(diagIndices)
    diagIndices <- repIndices == seqIndices
    rowIndices <- c(repIndices[ diagIndices],
                    repIndices[!diagIndices],
                    seqIndices[!diagIndices])
    colIndices <- c(seqIndices[ diagIndices],
                    seqIndices[!diagIndices],
                    repIndices[!diagIndices])
    vals <- numeric(length(rowIndices))
    vals <- c(diagVals, offDiagVals)
    valIndices <- c(1:matSize,
                    matSize + (1:sum(!diagIndices)),
                    matSize + (1:sum(!diagIndices)))
    ans <- repSparse(rowIndices, colIndices, valIndices, vals)
    class(ans) <- c("repSparseSymm", class(ans))
    return(ans)
}

setOldClass("repSparseSymm")
setIs("repSparseSymm", "repSparse")


##' Repeated sparse matrix with compound symmetry
##'
##' @param diagVal value for the diagonal
##' @param offDiagVal value for the off-diagonal
##' @param matSize size of the resulting matrix
##' @family repSparseSpecial
##' @export
##' @examples
##' (xCompSymm <- repSparseCompSymm(1, -0.2, 5))
repSparseCompSymm <- function(diagVal, offDiagVal, matSize) {
    ## MATNAME: Compound symmetry
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

setOldClass("repSparseCompSymm")
setIs("repSparseCompSymm", "repSparse")


##' Repeated sparse matrix with only one non-zero value off the
##' diagonal
##'
##' @param diagVal unique value for the diagonal
##' @param offDiagVal unique value for the off-diagonal
##' @param offDiagInds indices for the two correlated objects
##' @param matSize size of the matrix
##' @family repSparseSpecial
##' @export
##' @examples
##' (xOneOffDiag <- repSparseOneOffDiag(1, -0.2, c(5, 2), 5))
repSparseOneOffDiag <- function(diagVal, offDiagVal, offDiagInds, matSize) {
    ## MATNAME: One off diagonal
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

setOldClass("repSparseOneOffDiag")
setIs("repSparseOneOffDiag", "repSparse")


##' Repeated sparse Cholesky factor leading to constant variance
##'
##' @param sdVal standard deviation of crossproduct of the result
##' factor
##' @param offDiagVals off diagonal values of the Cholesky factor
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xConstVarChol <- repSparseConstVarChol(3, rnorm(choose(5, 2))))
repSparseConstVarChol <- function(sdVal, offDiagVals) {
    ## MATNAME: Constant variance Cholesky
    matSize <- nChoose2Inv(length(offDiagVals))
    diagIndices <- 1:matSize
    rowIndices <- rep(diagIndices, diagIndices)
    colIndices <- sequence(diagIndices)
    diagIndices <- rowIndices == colIndices
    vals <- numeric(length(diagIndices))
    innProd <- (sdVal^2) - c(0, tapply(offDiagVals^2, rowIndices[!diagIndices], sum))
    if(any(innProd < 0L)) stop("resulting matrix not positive definite")
    diagVals <- sqrt(innProd)
    ## a little DRY
    vals[ diagIndices] <- diagVals
    vals[!diagIndices] <- offDiagVals

    ans <- repSparse(rowIndices, colIndices,
                     seq_along(vals), vals,
                     mkConstVarCholTrans(c(sdVal, offDiagVals)))
    class(ans) <- c("repSparseConstVarChol", class(ans))
    return(ans)
}

setOldClass("repSparseConstVarChol")
setIs("repSparseConstVarChol", "repSparse")


##' Repeated sparse Cholesky factor of a correlation matrix
##'
##' @param offDiagPars parameters determining the off-diagonal of the
##' Cholesky factor
##' @family repSparseSpecial
##' @export
##' @examples
##' set.seed(1)
##' (xCorMatChol <- repSparseCorMatChol(rnorm(choose(5, 2))))
repSparseCorMatChol <- function(offDiagPars) {
    ## MATNAME: Correlation matrix Cholesky
    matSize <- nChoose2Inv(length(offDiagPars))
    diagIndices <- 1:matSize
    rowIndices <- rep(diagIndices, diagIndices)
    colIndices <- sequence(diagIndices)
    diagIndices <- rowIndices == colIndices
    vals <- numeric(length(diagIndices))
    splitPars <- split(offDiagPars, rowIndices[!diagIndices])
    offDiagValsList <- lapply(splitPars, function(xx) {
        xx <- xx^2
        sqrt(xx / (sum(xx) + 1))
    })
    innProd <- 1 - c(0, sapply(offDiagValsList, function(xx) sum(xx^2)))
    if(any(innProd < 0L)) stop("resulting matrix not positive definite") ## shouldn't ever happen
    diagVals <- sqrt(innProd)
    ## a little DRY
    vals[ diagIndices] <- diagVals
    vals[!diagIndices] <- sign(offDiagPars) * unlist(offDiagValsList)

    ans <- repSparse(rowIndices, colIndices,
                     seq_along(vals), vals,
                     mkCorMatCholTrans(offDiagPars))
    class(ans) <- c("repSparseCorMatChol", class(ans))
    return(ans)
}

setOldClass("repSparseCorMatChol")
setIs("repSparseCorMatChol", "repSparse")


##' Repeated sparse diagonal covariance matrix with a covariate
##' determining the diagonal
##'
##' @param varPars vector of variance parameters (one per level of
##' \code{grpFac})
##' @param covariate covariate
##' @param grpFac a grouping factor (or anything coercible to
##' \code{numeric} really) with the number of levels equal to the
##' length of \code{varPars}
##' @param mkTransFunc a function taking arguments \code{varPars},
##' \code{covariate}, and \code{grpFac} for constructing a
##' \code{trans} function (see \code{\link{mkVarExpTrans}} for an
##' example).
##' @family repSparseSpecial
##' @export
##' @examples
##' (xVarWithCovariate <- repSparseVarWithCovariate(rep(1, 5), 0.05*(4:0)))
repSparseVarWithCovariate <- function(varPars, covariate, grpFac,
                                      mkTransFunc = mkVarExpTrans) {
    ## MATNAME: Structured covariance matrix
    if(missing(grpFac)) grpFac <- factor(seq_along(covariate))
    if(missing(covariate)) covariate <- NA
    trans <- mkTransFunc(varPars, covariate, grpFac)
    matSize <- length(grpFac)
    vals <- trans(varPars)
    ans <- repSparse(1:matSize, 1:matSize, 1:matSize, vals,
                     trans = trans, Dim = c(matSize, matSize))
    class(ans) <- c("repSparseVarWithCovariate", class(ans))
    return(ans)
}

setOldClass("repSparseVarWithCovariate")
setIs("repSparseVarWithCovariate", "repSparse")



##' Cholesky factor of a repeated sparse covariance matrix obeying
##' exponential distance decay in covariance
##'
##' @param distObj distance matrix object
##' @param cutOffDist maximum distance with nonzero correlations
##' @family repSparseSpecial
##' @export
##' @examples
##' (xExpChol <- repSparseExpChol(dist(matrix(rnorm(10), 5, 2))))
repSparseExpChol <- function(distObj, cutOffDist = Inf) {
    ## MATNAME: Exponential distance Cholesky
    matSize <- as.integer(attr(distObj, "Size"))
    ri <- c(0:(matSize - 1L), matSize - rev(sequence(1:(matSize - 1))))
    ci <- c(1:matSize, rep(1:(matSize - 1), (matSize - 1):1)) - 1L
    va <- c(rep(0, matSize), as.numeric(distObj))

    ## FIXME: DRY -- order takes more than one argument?
    ord <- order(ri, decreasing = FALSE)
    ri <- ri[ord]; ci <- ci[ord]; va <- va[ord]
    ord <- order(ci, decreasing = FALSE)
    ri <- ri[ord]; ci <- ci[ord]; va <- va[ord]

    symmObj <- new("dsCMatrix", uplo = "L",
                   i = ri, p = ind2point(ci, matSize),
                   x = va,
                   Dim = rep(matSize, 2))

    vecDist <- symmObj@x
    symmObj@x <- exp(-vecDist)
    cholObj <- Cholesky(symmObj)

    cholMat <- as(cholObj, "sparseMatrix")
    ans <- repSparse(rowInds = cholMat@i + 1L,
                     colInds = point2ind(cholMat@p),
                     valInds = seq_along(vecDist),
                     vals = cholMat@x,
                     trans = mkExpCholTrans(1, cholObj, symmObj, vecDist))
    class(ans) <- c("repSparseExpChol", class(ans))
    return(ans)
}

setOldClass("repSparseExpChol")
setIs("repSparseExpChol", "repSparse")

##' Construct a repeated sparse upper Cholesky factor from an
##' \code{nlme}-style \code{corStruct} object
##'
##' @param object a \code{corStruct} object
##' @param sig initial standard deviation
##' @family repSparseSpecial
##' @export
##' @examples
##' (xCorFactor <- repSparseCorFactor(Initialize(corAR1(0.5, form = ~ 1 | Subject), data = Orthodont)))
repSparseCorFactor <- function(object, sig = 1) {
    ## MATNAME: Cholesky from corStruct object
    sigExists <- !is.null(sig)
    coefExists <- length(coef(object)) != 0
    
    corFac <- corFactor(object)
    lens <- Dim(object)$len
    vecLens <- 2 * choose(lens, 2) + lens
    vecList <- lme4ord:::subRagByLens(corFactor(object), vecLens)

    invList <- lme4ord:::mapplyInvList(vecList, lens)
    upperInds <- lapply(invList, upper.tri)

    oneBlock <- length(invList) == 1L
        ## ans <- repSparseTri(diagVals = diag(invList[[1]]),
        ##                     offDiagVals = invList[[1]][upperInds[[1]]],
        ##                     low = FALSE))

    ans <- .bind(mapply(repSparseTri,
                        lapply(invList, diag),
                        mapply("[", invList, upperInds, SIMPLIFY = FALSE),
                        MoreArgs = list(low = FALSE),
                        SIMPLIFY = FALSE), "diag")

    if(!coefExists) {
        ans <- resetTransConst(ans)
        if(!sigExists) {
            class(ans) <- c("repSparseCorFactor", class(ans))
            return(ans)
        }
    }
    if(coefExists) {
        diagIndices <- lapply(lens, seq, from = 1, by = 1)
        rowIndices <- mapply(rep, diagIndices, diagIndices, SIMPLIFY = FALSE)
        colIndices <- lapply(diagIndices, sequence)
        diagIndices <- mapply("==", rowIndices, colIndices, SIMPLIFY = FALSE)
        
        transEnv <- environment(ans$trans)
        list4env <- list(object = object,
                         init = coef(object),
                         lens = lens,
                         diagIndices = diagIndices,
                         rowIndices = rowIndices,
                         colIndices = colIndices,
                         vecLens = vecLens,
                         upperInds = upperInds)
        if(oneBlock) list4env$parList <- list(getInit(ans))
        list2env(list4env, transEnv)

        ans$trans <- local({
            function(matPars) {
                coef(object) <- matPars
                vecList <- subRagByLens(corFactor(object), vecLens)
                invList <- mapplyInvList(vecList, lens)
                diagVals <- lapply(invList, diag)
                upperVals <- mapply("[", invList, upperInds, SIMPLIFY = FALSE)
                for(i in seq_along(diagVals)) {
                    parList[[i]][ diagIndices[[i]]] <-  diagVals[[i]]
                    parList[[i]][!diagIndices[[i]]] <- upperVals[[i]]
                }
                unlist(parList)
            }
        }, transEnv)
    }

    if(sigExists) ans <- scalarMult(ans, sig)
    class(ans) <- c("repSparseCorFactor", class(ans))
    return(ans)
}

setOldClass("repSparseCorFactor")
setIs("repSparseCorFactor", "repSparse")

## ----------------------------------------------------------------------
## Random repeated sparse matrices
## ----------------------------------------------------------------------

##' Random repeated sparse matrices
##'
##' @param nrows,ncols numbers of rows and columns
##' @param nvals number of values
##' @param nnonzeros number of nonzero elements
##' @param rfunc random number function
##' @param ... dots
##' @export
rRepSparse <- function(nrows, ncols, nvals, nnonzeros, rfunc = rnorm, ...) {
    ## Random repeated sparse matrix
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
setMethod("chol", signature(x = "repSparse"), {
    function(x, ...) {
        ans <- as.repSparse(as(chol(as.matrix(x, sparse = TRUE)), "dgCMatrix"))
        class(ans) <- c("repSparseChol", class(x))
        return(ans)
    }
})

##' @rdname chol
##' @export
setMethod("chol", signature(x = "repSparseOneOffDiag"), {
    function(x, ...) {
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
        class(ans) <- c("repSparseChol", class(x))
        return(ans)
    }
})

##' @rdname chol
##' @export
##' @examples
##' x <- repSparseCompSymm(1.2, -0.11, 4)
##' as.matrix(chol(x))
##' t(chol(as.matrix(x)))
setMethod("chol", signature(x = "repSparseCompSymm"), {
    function(x, ...) {
        n <- nrow(x)
        trans <- mkCholCompSymmTrans(x$vals, n)
        sa <- seq_len(n)
        ii <- rep.int(1:(n-1), (n-1):1)
        ri <- c(sa, unlist(lapply(2:n, ":", n)))
        ci <- c(sa, ii)
        vi <- c(sa, ii + n)
        ans <- repSparse(ri, ci, vi,
                         trans(x$vals), trans = trans,
                         Dim = dim(x))
        class(ans) <- c("repSparseChol", class(x))
        return(ans)
    }
})


