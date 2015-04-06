##' Set model matrix slice and relative covariance factor block for a
##' random effects term
##' 
##' @param object a \code{\link{reTrmStruct}} object
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
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nlevels(object$grpFac)))),
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
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nlevels(object$grpFac)))),
                     class = class(object)))
}
