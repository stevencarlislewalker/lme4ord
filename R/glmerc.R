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
                parInds = parInds, X = parsedForm$X, y = parsedForm$y,
                lower = lower)
    class(ans) <- "glmerc"
    return(ans)
}

##' @export
fixef.glmerc <- function(object, ...) .fixef(object$opt$par, object$parInds)

##' @export
covar.glmerc <- function(object, ...) .covar(object$opt$par, object$parInds)

##' @export
loads.glmerc <- function(object, ...) stop("covariance over levels models do not have loadings")

##' @export
covarByTerms <- function(object, ...) {
    nCovar <- object$parsedForm$nCovar
    if(length(nCovar) == 1L) return(list(covar(object)))
    inds <- covarInds(object$parsedForm$nCovar)
    lapply(inds, function(i) covar(object)[i])
}

##' @export
covarInds <- function(nCovar) {
    ans <- lapply(nCovar, seq, from = 1, by = 1)
    for(i in 2:length(nCovar)) ans[[i]] <- ans[[i]] + max(ans[[i-1]])
    return(setNames(ans, names(nCovar)))
}

##' @export
VarCorr.glmerc <- function(x, ...) {
    tmm <- getMEc(x, "TmodMat")
    cnms <- getMEc(x, "cnms")
    lapply(lapply(tmm, as.matrix), crossprod)
}

##' @export
vcov.glmerc <- function(object, justFixef = TRUE, ...) {
    optPar <- object$opt$par
    ans <- solve(0.5 * lme4:::deriv12(object$dfun,
                                      optPar)$Hessian)
    dimnames(ans) <- rep(list(names(optPar)), 2)
    if(justFixef) {
        dims <- object$parInds$fixef
        ans <- ans[dims, dims]
    }
    return(ans)
}

##' @export
print.glmerc <- function(x, ...) {
    cat("\nGeneralized linear mixed model\nwith covariance amongst grouping factor levels\n")
    cat("----------------------------------------------\n\n")

    cat("Fixed effects\n")
    cat("-------------\n\n")
    print(cbind(Estimate = fixef(modNest),
                `Std. Error` = sqrt(diag(vcov(modNest)))))

    cat("\n\nRandom effects (co)variance\n")
    cat("---------------------------\n\n")
    print(VarCorr(x))
}
