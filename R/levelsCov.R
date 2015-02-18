##' Parse mixed model formula for levels-covariance model
##'
##' @param formula mixed model formula
##' @param data data
##' @param family family
##' @param covList list of covariance matrices for random effects
##' grouping factors in \code{formula}
##' @param ... not used yet
##' @export
levelsCovFormula <- function(formula, data = NULL, family = binomial, covList, ...) {
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")
    for(i in seq_along(reTrmsList)) {
        reTrmsList[[i]]$covMat <- covList[[reTrmsList[[i]]$grpName]]
        if(is.null(reTrmsList[[i]]$covMat)) {
            nli <- nlevels(reTrmsList[[i]]$grpFac)
            reTrmsList[[i]]$covMat <- diag(1, nli, nli)
        }
    }
    modMats <- lapply(reTrmsList, "[[", "modMat")
    grpFacs <- lapply(reTrmsList, "[[", "grpFac")
    grpFacConst <- rep(list(as.factor(rep("const", nrow(data)))), length(grpFacs))
    covMats <- lapply(reTrmsList, "[[", "covMat")
    covMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))
    re <- mkTemplateReTrms(modMats,
                           grpFacs, grpFacConst,
                           covMats, covMatConst)
    return(c(list(X = model.matrix(nobars(formula), data),
                  y = model.response(model.frame(nobars(formula), data))), re))
    ## initPars <- list(covar = re$theta, fixef = rep(0, ncol(X)), loads = NULL)
    ## parInds <- lapply(initPars, seq_along)
    ## mkGeneralGlmerDevfun(y, X, re$Zt, re$Lambdat,
    ##                      rep(1, length(y)), rep(0, length(y)),
    ##                      initPars, parInds,
    ##                      re$mapToCovFact, function(loads) NULL)
}

##' Get model matrix and grouping factor
##' 
##' @param bar random effect language object
##' @param fr model frame
##' @return list with model matrix and grouping factor
##' @export
getModMatAndGrpFac <- function(bar, fr) {
    ## based on mkBlist
    
    fr <- factorize(bar, fr)
    grpLang <- bar[[3]]
    linFormLang <- bar[[2]]
    nm <- deparse(grpLang)
    ## try to evaluate grouping factor within model frame ...
    if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac),
                                               list(fac = grpLang)), fr),
                error=function(e) NULL)))
        stop("couldn't evaluate grouping factor ",
             nm," within model frame:",
             " try adding grouping factor to data ",
             "frame explicitly if possible",call.=FALSE)
    if (all(is.na(ff)))
        stop("Invalid grouping factor specification, ",
             nm, call. = FALSE)
    mm <- model.matrix(eval(substitute( ~ foo, list(foo = linFormLang))), fr)
    return(list(modMat = mm, grpFac = ff, grpName = nm))
}
