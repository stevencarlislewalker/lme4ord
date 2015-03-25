generalGlmer <- function(formula, data, ...) {
    generalParseFormula(formula, data, ...)
    
}

##' Parse a mixed model formula
##'
##' @param formula mixed model formula
##' @param data an object coercible to data frame
##' @param ... additional parameters to \code{\link{as.data.frame}}
##' @export
generalParseFormula <- function(formula, data, ...) {
    data <- as.data.frame(data, ...)
                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

    return(list(response = model.response(model.frame(nobars(formula), data)),
                fixed    = model.matrix(nobars(formula), data),
                random   = reTrmsList))
}

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
