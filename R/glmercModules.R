##' Parse mixed model formula for levels-covariance model
##'
##' @param formula mixed model formula
##' @param data data
##' @param family family
##' @param covList list of covariance matrices for random effects
##' grouping factors in \code{formula}
##' @param ... not used yet
##' @export
glmercFormula <- function(formula, data = NULL, family = binomial, covList, ...) {

                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

                                        # append the covariance
                                        # matrices over the levels
                                        # associated with each random
                                        # effects term
    for(i in seq_along(reTrmsList)) {
        reTrmsList[[i]]$covMat <- covList[[reTrmsList[[i]]$grpName]]
        if(is.null(reTrmsList[[i]]$covMat)) { # if no cov mat, use
                                              # identity
            nli <- nlevels(reTrmsList[[i]]$grpFac)
            reTrmsList[[i]]$covMat <- diag(1, nli, nli)
        }
    }

                                        # extract the model matrices,
                                        # grouping factors, and cov
                                        # mats
    modMats <- lapply(reTrmsList, "[[", "modMat")
    grpFacs <- lapply(reTrmsList, "[[", "grpFac")
    covMats <- lapply(reTrmsList, "[[", "covMat")

                                        # set up constant grouping
                                        # factors and 1 by 1 cov mats
                                        # for mkTemplateReTrms (this
                                        # is where we could add
                                        # nesting)
    grpFacConst <- rep(list(as.factor(rep("const", nrow(data)))), length(grpFacs))
    covMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))

                                        # construct the full random
                                        # effects model matrix,
                                        # relative covariance factor,
                                        # etc.
    re <- mkTemplateReTrms(modMats,
                           grpFacs, grpFacConst,
                           covMats, covMatConst)

                                        # organize output value
    return(c(re,
             list(X = model.matrix(nobars(formula), data),
                  y = model.response(model.frame(nobars(formula), data)),
                  flist = simplifyFacList(grpFacs),
                  modMats = modMats)))
}

##' @export
getMEc <- function(object, name = c("X", "Z", "Zt", "y", "TmodMat", "cnms")) {
    if (missing(name)) 
        stop("'name' must not be missing")
    if (length(name <- as.character(name)) > 1) {
        names(name) <- name
        return(lapply(name, getMEc, object = object))
    }
    cnms <- lapply(object$parsedForm$modMats, colnames)
    TmodMat <- object$parsedForm$TmodMat
    
    for(i in 1:length(TmodMat)) {
        TmodMat[[i]]@x <- covarByTerms(object)[[i]]
        dimnames(TmodMat[[i]]) <- rep(cnms[i], 2)
    }
    switch(name,
           X = object$X,
           Z = t(object$parsedForm$Zt),
           Z = object$parsedForm$Zt,
           y = object$y,
           TmodMat = TmodMat,
           cnms = cnms)
}

##' Make function for computing deviance profiles
##' 
##' @param mod \code{\link{glmerc}} object
##' @param whichPar which parameter to profile
##' @return function for computing the profiled deviance over a single
##' parameter
##' @export
mkPfun <- function(mod, whichPar) {
    local({
        dfun <- mkNfun(mod)
        optPar <- mod$opt$par
        whichPar <- whichPar
        op <- optPar[-whichPar]
        function(par) {
            fn <- function(op) {
                parNow <- numeric(length(op) + 1)
                parNow[whichPar] <- par
                parNow[-whichPar] <- op
                dfun(parNow)
            }
            opt <- bobyqa(op, fn, lower = mod$lower[-whichPar],
                          control = list(iprint = 0L, maxfun = 500))
            op <<- opt$par # FIXME: assign to local env
            return(opt$fval)
        }
    })
}

##' Make function for computing deviance slices
##' 
##' @param mod \code{\link{glmerc}} object
##' @param whichPar which parameter to slice
##' @return function for computing 1D slices through the deviance at
##' the optimum
##' @export
mkSfun <- function(mod, whichPar) {
    local({
        dfun <- mkNfun(mod)
        whichPar <- whichPar
        par <- mod$opt$par
        function(par1D) {
            par[whichPar] <- par1D
            return(dfun(par))
        }
    })
}

##' Make function for computing normalized deviance
##'
##' @param mod
##' @return function for computing the normalized deviance, which
##' equals zero at the optimum
##' @export
mkNfun <- function(mod) {
    local({
        dfun <- mod$dfun
        optFunVal <- dfun(mod$opt$par)
        function(par) dfun(par) - optFunVal
    })
}
