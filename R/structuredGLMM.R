## ----------------------------------------------------------------------
## Structured GLMM
##
## uses:  strucParseFormula.R, mkGeneralGlmerDevfun.R, outputModule.R
## ----------------------------------------------------------------------

##' Generalized linear mixed model with structured (co)variance terms
##'
##' Fits a generalized linear mixed model with structured (co)variance
##' terms by the following fully modularized steps:
##' \describe{
##' \item{\code{\link{strucParseFormula}}}{parses the mixed model
##' formula with structured terms}
##' \item{\code{\link{mkGeneralGlmerDevfun}}}{constructs a generalized
##' linear mixed model deviance function}
##' \item{\code{bobyqa}}{optimizer the deviance function}
##' \item{\code{\link{mkStrucGlmer}}}{constructs an object of
##' \code{\link{strucGlmer-class}}}}
##'
##' @param formula extended mixed model formula
##' @param data,family,weights,offset,etastart,devfunOnly see \code{\link{glmer}} 
##' @param addArgs list of additional arguments to pass to
##' @param optVerb verbose
##' @param optMaxit maximum number of iterations for the nonlinear
##' optimizers
##' @param ... further arguments to \code{\link{mkGeneralGlmerDevfun}}
##' @export
strucGlmer <- function(formula, data, family, addArgs = list(),
                       optVerb = 0L, optMaxit = 10000,
                       weights = NULL, offset = NULL, etastart = NULL,
                       devfunOnly = FALSE,
                         ...) {

    mc <- match.call()
    
    cat("\nConstructing vectors and matrices...\n")
    parsedForm <- strucParseFormula(formula, data, addArgs, ...)

    cat("\nConstructing deviance function...\n")
    dfun <- mkGeneralGlmerDevfun(y = parsedForm$response,
                                 X = parsedForm$fixed,
                                 Zt      = as(parsedForm$Zt,      "dgCMatrix"),
                                 Lambdat = as(parsedForm$Lambdat, "dgCMatrix"),
                                 weights   = weights,
                                 offset    = offset,
                                 etastart  = etastart,
                                 initPars  = parsedForm$initPars,
                                 parInds   = parsedForm$parInds,
                                 mapToCovFact = parsedForm$mapToCovFact,
                                 mapToModMat  = parsedForm$mapToModMat,
                                 devfunEnv = parsedForm$devfunEnv,
                                 family = family,
                                 Lind = parsedForm$Lambdat$valInds,
                                 ...)
    if(devfunOnly) return(dfun)
    
    cat("\nOptimizing deviance function...\n")
    opt <- minqa::bobyqa(parsedForm$initPars, dfun,
                         lower = parsedForm$lower,
                         upper = parsedForm$upper,
                         control =
                         list(iprint = optVerb,
                              maxfun = optMaxit,
                              rhobeg = 0.0002,
                              rhoend = 2e-7))
    
    cat("\nPreparing output...\n")
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}

