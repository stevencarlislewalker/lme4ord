##' Generalized linear mixed model with custom covariance over grouping factor levels
##'
##' When fitting generalized linear mixed model with
##' \code{\link{glmer}}, the covariance structure over the grouping
##' factor levels is assumed to be an identity matrix.  That is,
##' random effects are sampled indepedently over the grouping factor
##' levels.  With \code{glmerc}, one may specify the covariance over
##' the levels, up to a fitted parameter.
##'
##' @param formula formula
##' @param data data
##' @param family family
##' @param covList covList
##' @param optControl optControl
##' @param ... ...
##' @export
glmerc <- function(formula, data = NULL, family = binomial, covList = list(),
                   optControl = list(iprint = 0L), ...) {

                                        # parse formula
    data <- as.data.frame(data)
    parsedForm <- glmercFormula(formula, data, covList = covList)

                                        # organize initial values
    covar <- parsedForm$covar
    fixef <- rep(0, ncol(parsedForm$X))
    initPars <- c(covar = covar,
                  fixef = fixef)
    parInds <- list(covar = seq_along(covar),
                    fixef = seq_along(fixef) + length(covar),
                    loads = NULL)

                                        # construct deviance function
    dfun <- mkGeneralGlmerDevfun(parsedForm$y, parsedForm$X,
                                 parsedForm$Zt, parsedForm$Lambdat,
                                 rep(1, nrow(data)), rep(0, nrow(data)),
                                 initPars, parInds,
                                 parsedForm$mapToCovFact, function(loads) NULL)

                                        # optimize deviance function
    dfun(initPars)
    lower <- ifelse(initPars, 0, -Inf)
    opt <- bobyqa(initPars, dfun, lower = lower,
                  control = list(iprint = 4L))
    names(opt$par) <- names(initPars)

                                        # organize return value
    ans <- list(opt = opt, parsedForm = parsedForm, dfun = dfun,
                parInds = parInds, X = X, y = y)
    class(ans) <- "glmerc"
    return(ans)
}

##' @export
fixef.glmerc <- function(object, ...) ans$opt$par[ans$parInds$fixef]
