##' Phylogenetic generalized linear mixed model
##'
##' @param formula formula
##' @param data data
##' @param family family
##' @param covList covList
##' @param optControl optControl
##' @param ... ...
##' @export
pglmer <- function(formula, data = NULL, family = binomial, covList = list(),
                   optControl = list(iprint = 0L), ...) {

    data <- as.data.frame(data)
    parsedForm <- levelsCovFormula(formula, data, covList = covList)
    X <- model.matrix(nobars(formula), df)
    fr <- model.frame(formula, df)
    y <- model.response(fr)
    covar <- parsedForm$covar
    fixef <- rep(0, ncol(X))
    initPars <- c(covar = covar,
                  fixef = fixef)
    parInds <- list(covar = seq_along(covar),
                    fixef = length(covar),
                    loads = NULL)
    dfun <- mkGeneralGlmerDevfun(y, X,
                                 parsedForm$Zt, parsedForm$Lambdat,
                                 rep(1, nrow(data)), rep(0, nrow(data)),
                                 initPars, parInds,
                                 parsedForm$mapToCovFact, function(loads) NULL,
                                 family = family())
    dfun(initPars)
    lower <- ifelse(covar, 0, -Inf)
    opt <- bobyqa(initPars, dfun, lower = lower,
                  control = list(iprint = 4L))
    names(opt$par) <- names(initPars)
    return(list(opt = opt, parsedForm = parsedForm, dfun = dfun,
                parInds = parInds, X = X, y = y))
}
