##' Optimize a structured GLMM deviance function
##'
##' @param optimizer nonlinear optimization function (currently only
##' minqa::bobyqa is allowed)
##' @param initPars initial parameter values
##' @param dfun structured GLMM deviance function (results of
##' \code{\link{strucMkDevfun}} or \code{\link{mkGeneralGlmerDevfun}})
##' @param lower,upper bounds on the parameter vector
##' @param optVerb,optMaxit controls for the optimizer (see
##' \code{\link{strucGlmer}})
##' @param ... not used
##' @export
strucOpt <- function(optimizer = minqa::bobyqa,
                     initPars, dfun,
                     lower, upper,
                     optVerb, optMaxit, ...) {

    if(!identical(optimizer, minqa::bobyqa))
        stop("minqa::bobyqa is currently the only available optimizer")

    try(optimizer(initPars, dfun, lower = lower, upper = upper,
                  control =
                  list(iprint = optVerb,
                       maxfun = optMaxit,
                       rhobeg = 0.0002,
                       rhoend = 2e-7)))
}

