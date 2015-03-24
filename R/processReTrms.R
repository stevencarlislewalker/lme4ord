##' Make general random effects terms
##'
##' @param reStructList list of random effects structure functions
##' @param parsedForm list parsed random effects terms
##' @param ... potential extra parameters
##' @rdname mkReTrms
##' @export
mkGeneralReTrms <- function(reStructList, parsedForm, ...) {
    trmList <- mapply(do.call, reStructList, parsedForm$random, SIMPLIFY = FALSE)
    trmList <- listTranspose(trmList)
         ZtBind <- as.repSparse(sort(.bind(trmList$Zt,      "row" )))
    LambdatBind <- as.repSparse(sort(.bind(trmList$Lambdat, "diag")))
    return(list(Zt = as(ZtBind, "dgCMatrix"),
                Lambdat = as(LambdatBind, "dgCMatrix"),
                ZtTrans = mkSparseTrans(ZtBind),
                LambdatTrans = mkSparseTrans(LambdatBind)))
}

##' Make transformation function for column-compressed sparse matrix
##'
##' @param object repeated sparse matrix object
##' @rdname mkReTrms
##' @export
mkSparseTrans <- function(object) {
    local({
        trans <- object$trans
        inds <- object$valInds
        function(matPars) trans(matPars)[inds]
    })
}

##' @param modMat raw random effects model matrix
##' @param grpFac grouping factor
##' @param grpName grouping factor name
##' @rdname mkReTrms
##' @export
mkIdentityReTrms <- function(modMat, grpFac, grpName) {
    Zt <- kr(t(as.repSparse(modMat)), as.repSparse(grpFac))
    return(list(Zt = resetTransConst(Zt),
                Lambdat = repSparseIdent(nlevels(grpFac))))
}
