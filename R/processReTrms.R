##' Get random effects term
##' 
##' @param modMat raw random effects model matrix
##' @param grpFac grouping factor
##' @param grpName grouping factor name
##' @rdname getReTrm
##' @export
getReTrm <- function(modMat, grpFac, grpName) {
    UseMethod("getReTrm")
}

##' @rdname getReTrm
##' @export
getReTrm.identity <- function(modMat, grpFac, grpName) {
    Zt <- kr(t(as.repSparse(modMat)), as.repSparse(grpFac))
    return(list(Zt = resetTransConst(Zt),
                Lambdat = repSparseIdent(nlevels(grpFac))))
}

