##' Make general random effects terms
##'
##' @param trmList list of terms
##' @param ... potential extra parameters
##' @rdname mkReTrms
##' @export
mkGeneralReTrms <- function(trmList, ...) {
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

##' @param trm random effects term
##' @rdname mkReTrms
##' @export
mkIdentityReTrms <- function(trm, ...) {
    Zt <- kr(t(as.repSparse(trm$modMat)), as.repSparse(trm$grpFac))
    return(list(Zt = resetTransConst(Zt),
                Lambdat = repSparseIdent(nlevels(trm$grpFac))))
}
