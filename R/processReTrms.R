##' Set model matrix slice and relative covariance factor block for a
##' random effects term
##' 
##' @param object a \code{reTrmStruct} object
##' @param addArgs named list of additional arguments
##' @rdname setReTrm
##' @export
setReTrm <- function(object, addArgs) {
    UseMethod("setReTrm")
}

##' @rdname setReTrm
##' @export
setReTrm.default <- function(object, addArgs) {
    Zt <- kr(t(as.repSparse(object$modMat)), as.repSparse(object$grpFac))
    nc <- ncol(object$modMat)
    nl <- nlevels(object$grpFac)
    template <- repSparseTri(rep(1, nc), rep(0, choose(nc, 2)), FALSE)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = rep(template, nl, type = "diag"))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.identity <- function(object, addArgs) {
    Zt <- kr(t(as.repSparse(object$modMat)), as.repSparse(object$grpFac))
    nl <- nlevels(object$grpFac)
    nc <- ncol(object$modMat)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nl * nc))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.edge <- function(object, addArgs) {
    addArgs <- eval(object$addArgs, addArgs)
    Jedge <- as(edgeTipIndicator(addArgs$phy), "sparseMatrix")
    Jspp <- as(object$grpFac, "sparseMatrix")
    Jt <- as.repSparse(Jedge %*% Jspp)
    Zt <- kr(t(as.repSparse(object$modMat)), Jt)
    nl <- nrow(Jedge)
    nc <- ncol(object$modMat)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nl * nc))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.cooccur <- function(object, addArgs) {
    addArgs <- eval(object$addArgs, addArgs)
    Jt <- as.repSparse(as(object$grpFac, "sparseMatrix"))
    Zt <- kr(t(as.repSparse(object$modMat)), Jt)
    nCovPars <- choose(ncol(object$modMat), 2)
    Tt <- t(repSparseCorMatChol(rep(0, nCovPars)))
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = rep(Tt, nrow(Jt), type = "diag")),
                       class = class(object))))
}

##' Simulate additional arguments
##'
##' @param object a \code{reTrmStruct} object
##' @param ... dots
##' @rdname simAddArgs
##' @export
simAddArgs <- function(object, ...) {
    UseMethod("simAddArgs")
}

##' @rdname simAddArgs
##' @export
simAddArgsList <- function(object, ...) {
    l... <- list(...)
    unlist(do.call(lapply, c(list(unname(object), simAddArgs), l...)), FALSE)
}

##' @rdname simAddArgs
##' @export
simAddArgs.default <- function(object, ...) list()

##' @param rtreeArgs arguments for \code{\link{rtree}}
##' @param compute.brlenArgs arguments for \code{\link{compute.brlen}}
##' @rdname simAddArgs
##' @export
simAddArgs.edge <- function(object, rtreeArgs = list(),
                            compute.brlenArgs = list(), ...) {
    namePhy <- as.character(object$addArgs$phy)
    phy <- do.call(rtree, c(list(nlevels(object$grpFac)), rtreeArgs))
    phy <- do.call(compute.brlen, c(list(phy), compute.brlenArgs))
    phy$tip.label <- unique(as.character(object$grpFac))
    return(setNames(list(phy), namePhy))
}


