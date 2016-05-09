## ----------------------------------------------------------------------
## Structured GLMM
##
## uses:  formulaParsing.R, devianceFunctionConstruction.R, outputAndInference.R
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
##' @param formula extended mixed model formula or
##' \code{\link{strucParseFormula}} object
##' @param data,family,weights,offset,etastart,devfunOnly see \code{\link{glmer}} 
##' @param addArgs list of additional arguments to pass to
##' @param modularVerb should the modular steps be announced during the fitting process?
##' @param optVerb verbose
##' @param optMaxit maximum number of iterations for the nonlinear
##' optimizers
##' @param numModularSteps number of modular steps to perform
##' (defaults to \code{4}, which is all of them). \code{1} returns the
##' \code{\link{strucParseFormula}} object, \code{2} the deviance
##' function, \code{3} the optimizer output, and \code{4} a
##' \code{strucGlmer} object.
##' @param parList potential named list of initial parameters
##' @param parsedForm \code{\link{strucParseFormula}} object
##' @param initPars initial parameter values (required for
##' \code{strucGlmer.function} method if missing \code{parsedForm})
##' @param lower,upper lower and upper bounds of parameter values
##' (only useful for \code{strucGlmer.function} method)
##' @param ... further arguments to \code{\link{mkGeneralGlmerDevfun}}
##' @importFrom minqa bobyqa
##' @export
##' @examples
##' cbpp$incidenceBySize <- with(cbpp, incidence/size)
##' gm <- strucGlmer(incidenceBySize ~ factAnal(0 + herd | period, nAxes = 1),
##'                  data = cbpp, family = binomial, weights = cbpp$size,
##'                  penLoads = mkPenLpNorm())
strucGlmer <- function(formula, ...) {
    UseMethod("strucGlmer")
}

##' @rdname strucGlmer
##' @export
strucGlmer.default <- function(formula, ...) {
    stop("no strucGlmer method for objects of class", class(formula))
}

##' @rdname strucGlmer
##' @export
strucGlmer.formula <- function(formula, data, family, addArgs = list(),
                               modularVerb = FALSE, optVerb = 0L, optMaxit = 10000,
                               weights = NULL, offset = NULL, etastart = NULL,
                               devfunOnly = FALSE, numModularSteps = 4L, parList = NULL,
                               ...) {

    mc <- match.call()
    numModularSteps <- round(numModularSteps)
    
    if(modularVerb) cat("\nConstructing vectors and matrices...\n")
    parsedForm <- strucParseFormula(formula, data, addArgs, parList = parList, ...)
    if(numModularSteps == 1L) return(parsedForm)

    if(modularVerb) cat("\nConstructing deviance function...\n")
    dfun <- strucMkDevfun(parsedForm, family = family, weights = weights,
                          offset = offset, etastart = etastart, ...)
    if(devfunOnly || (numModularSteps == 2L)) return(dfun)
    
    if(modularVerb) cat("\nOptimizing deviance function...\n")
    opt <- strucOpt(minqa::bobyqa, parsedForm$initPars, dfun,
                    parsedForm$lower, parsedForm$upper, optVerb, optMaxit)
    if(numModularSteps == 3L) return(opt)

    if(modularVerb) cat("\nPreparing output...\n")
    out <- mkStrucGlmer(opt, parsedForm, dfun, mc)
    if(numModularSteps == 4L) return(out)

    stop(paste("cannot run", numModularSteps, "modular steps"))
}

##' @rdname strucGlmer
##' @export
strucGlmer.strucParseFormula <- function(formula, family, optVerb = 0L, optMaxit = 10000,
                                         weights = NULL, offset = NULL, etastart = NULL,
                                         ...) {
    mc <- match.call()
    parsedForm <- formula
    if(is.null(fam <- stats::family(parsedForm))) fam <- family
    dfun <- strucMkDevfun(parsedForm, family = fam, weights = weights,
                          offset = offset, etastart = etastart, ...)

    opt <- strucOpt(minqa::bobyqa, parsedForm$initPars, dfun,
                    parsedForm$lower, parsedForm$upper, optVerb, optMaxit)
    
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}

##' @rdname strucGlmer
##' @export
strucGlmer.function <- function(formula, optVerb = 0L, optMaxit = 10000,
                                parsedForm, initPars, lower, upper, ...) {
    mc <- match.call()
    if(!missing(parsedForm)) {
        initPars <- parsedForm$initPars
        lower <- parsedForm$lower
        upper <- parsedForm$upper
    } else {
        if(missing(initPars)) stop("need initPars if no parsedForm")
        if(missing(lower))    stop("need lower if no parsedForm")
        if(missing(upper))    stop("need upper if no parsedForm")
    }
    dfun <- formula
    opt <- strucOpt(minqa::bobyqa, initPars, dfun,
                    lower, upper, optVerb, optMaxit)
    if(missing(parsedForm)) return(opt)
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

    dfun <- strucMkDevfun(parsedForm, family = family(object), weights = weights,
                          offset = offset, etastart = etastart, ...)

    opt <- strucOpt(minqa::bobyqa, parsedForm$initPars, dfun,
                    parsedForm$lower, parsedForm$upper, optVerb, optMaxit)
    
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}

