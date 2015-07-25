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
##' @param parList potential named list of initial parameters
##' @param ... further arguments to \code{\link{mkGeneralGlmerDevfun}}
##' @importFrom minqa bobyqa
##' @export
##' @examples
##' cbpp$incidenceBySize <- with(cbpp, incidence/size)
##' gm <- strucGlmer(incidenceBySize ~ factAnal(0 + herd | period, nAxes = 1),
##'                  data = cbpp, family = binomial, weights = cbpp$size,
##'                  penLoads = mkPenLpNorm())
strucGlmer <- function(formula, data, family, addArgs = list(),
                       optVerb = 0L, optMaxit = 10000,
                       weights = NULL, offset = NULL, etastart = NULL,
                       devfunOnly = FALSE, parList = NULL,
                         ...) {

    mc <- match.call()
    
    cat("\nConstructing vectors and matrices...\n")
    parsedForm <- strucParseFormula(formula, data, addArgs, parList = parList, ...)

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
    opt <- try(minqa::bobyqa(parsedForm$initPars, dfun,
                             lower = parsedForm$lower,
                             upper = parsedForm$upper,
                             control =
                             list(iprint = optVerb,
                                  maxfun = optMaxit,
                                  rhobeg = 0.0002,
                                  rhoend = 2e-7)))
    
    cat("\nPreparing output...\n")
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}

##' @param object \code{strucGlmer} object
##' @method update strucGlmer
##' @rdname strucGlmer
##' @export
update.strucGlmer <- function(object, formula, parList,
                              data = NULL, addArgs = list(),
                              optVerb = 0L, optMaxit = 10000,
                              weights = NULL, offset = NULL,
                              etastart = NULL, ...) {
    mc <- match.call()
    
    if(!missing(formula)) {
        if(is.null(data)) stop("data not saved in strucGlmer objects,\n",
                               "therefore data must be provided when updating the formula")
        parsedForm <- strucParseFormula(update(formula(object), formula), data, addArgs, ...)
    } else {
        parsedForm <- object$parsedForm
    }
    if(!missing(parList)) {
        parsedForm <- update(parsedForm, parList)
    }
    if(is.null(weights)) weights <- weights(object)
    if(is.null(offset)) offset <- getOffset(object)
    
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
                                 family = family(object),
                                 Lind = parsedForm$Lambdat$valInds,
                                 ...)

    opt <- try(minqa::bobyqa(parsedForm$initPars, dfun,
                             lower = parsedForm$lower,
                             upper = parsedForm$upper,
                             control =
                             list(iprint = optVerb,
                                  maxfun = optMaxit,
                                  rhobeg = 0.0002,
                                  rhoend = 2e-7)))
    
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}



