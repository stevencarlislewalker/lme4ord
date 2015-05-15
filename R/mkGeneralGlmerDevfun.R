##' Make a general (Laplace) glmer deviance function
##'
##' @param y response vector
##' @param X fixed effects model matrix
##' @param Zt transposed random effects model matrix
##' @param Lambdat relative covariance factor
##' @param weights prior weights
##' @param offset known additive offsets
##' @param initPars initial values for the parameter vector
##' @param parInds a named list of vectors for identifying each of the
##' three types of parameters in \code{initPars}, one vector for each
##' of the three types of parameters: \code{list(covar = ., fixef = .,
##' loads = .)} where the \code{.}'s are vectors of indices
##' @param mapToCovFact function taking the \code{covar} parameters
##' and returning the values of the non-zero elements of
##' \code{Lambdat} (i.e. the \code{x} slot of \code{Lambdat})
##' @param mapToModMat function taking the \code{loads} parameters and
##' returning the values of the non-zero elements of \code{Zt}
##' (i.e. the \code{x} slot of \code{Zt})
##' @param family \code{\link{family}}
##' @param tolPwrss tolerance for penalized weighted residual sum of
##' squares
##' @param verbose verbose
##' @param pureR should the PIRLS algorithm be run in
##' \code{lme4pureR}?
##' @param Lind optional indices mapping covariance parameters to the
##' non-zero elements of the relative covariance factor
##' @return a deviance function with an environment
##' @importFrom lme4pureR pirls
##' @importFrom lme4 GHrule
##' @export
mkGeneralGlmerDevfun <- function(y, X, Zt, Lambdat,
                                 weights = NULL, offset = NULL,
                                 etastart = NULL, mustart = NULL,
                                 devfunEnv,
                                 initPars, parInds,
                                 mapToCovFact = NULL,
                                 mapToModMat = NULL,
                                 mapToWeights = NULL,
                                 family,
                                 tolPwrss = 1e-6,
                                 verbose = 0L, pureR = FALSE,
                                 Lind = NULL) {
    
    if(pureR) {
        stop("pure R implementation not currently working")
        return(pirls(X, y, Zt, Lambdat, mapToCovFact,
                     initPars, weights, offset,
                     family = family, tol = tolPwrss))
    }
    if(is.null(Lind)) {
        theta <- Lambdat@x
        Lind <- seq_along(Lambdat@x)
    } else {
        theta <- initPars[parInds$covar]
    }

    family <- fixFamily(family)

    initializeEnv <- initializeResp(y, etastart, mustart, offset, weights, family)

                                        # handle potential missing
                                        # arguments
    if(missing(parInds)) {
        if(is.recursive(initPars)) {
            parInds <- mkParInds(initPars)
        } else if(is.relistable(initPars)) {
            parInds <- mkParInds(relist(initPars))
        } else {
            stop("can't make parInds, please supply it")
        }
    }
    initPars <- unlist(initPars)

    if(missing(devfunEnv)) devfunEnv <- new.env()
    if(is.null(weights))   weights   <- rep(1, length(y))
    if(is.null(offset))    offset    <- rep(0, length(y))
    if(is.null(etastart))  etastart <- family$linkfun(initializeEnv$mustart)

    devfunList <- list(Lind = Lind,
                       pp = lme4:::merPredD$new(
                           X = X, Zt = Zt,
                           Lambdat = Lambdat,
                           Lind = Lind,
                           theta = as.double(theta),
                           n = nrow(X)),
                       resp = lme4:::glmResp$new(
                           y = y, family = family,
                           weights = weights),
                       lp0 = etastart,
                       baseOffset = offset,
                       tolPwrss = tolPwrss,
                       GQmat = GHrule(1), ## always Laplace approx
                       compDev = TRUE,
                       fac = NULL,
                       verbose = verbose,
                       setCovar = !is.null(parInds$covar),
                       setLoads = !is.null(parInds$loads),
                       setWeigh = !is.null(parInds$weigh),
                       setFixef = !is.null(parInds$fixef),
                       mapToCovFact = mapToCovFact,
                       mapToModMat  = mapToModMat,
                       mapToModMat  = mapToWeights,
                       parInds = parInds)

    devfun <- function(pars) {
        resp$setOffset(baseOffset)
        resp$updateMu(lp0)
        if(setCovar) pp$setTheta(as.double(mapToCovFact(pars[parInds$covar])))
        if(setLoads) pp$setZt(as.double(mapToModMat(pars[parInds$loads])))
        if(setWeigh) resp$setWeights(as.double(mapToWeights(pars[parInds$weigh])))
        spars <- as.numeric(pars[parInds$fixef])
        offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
        resp$setOffset(offset)
        p <- lme4:::glmerPwrssUpdate(pp, resp, tolPwrss, GQmat,
                                     compDev, fac, verbose)
        resp$updateWts()
        p
    }
    
    environment(devfun) <- list2env(devfunList, envir = devfunEnv)

                                        # initialize weights etc ...
    devfunEnv$resp$updateMu(etastart)
    devfunEnv$resp$updateWts()
    devfun(initPars)

    return(devfun)
}

##' Get parameters
##'
##' @param object \code{lme4ord} fitted model object
##' @rdname pars
##' @export
covar <- function(object, ...) UseMethod("covar")

##' @rdname pars
##' @export
loads <- function(object, ...) loadings(object)

##' @rdname pars
##' @export
pars <- function(object, ...) UseMethod("pars")

.covar <- function(pars, ind) pars[ind$covar]
.fixef <- function(pars, ind) pars[ind$fixef]
.loads <- function(pars, ind) pars[ind$loads]

##' @param parList named list of parameters with possible names:
##' (\code{covar}, \code{fixef}, \code{weigh}, \code{loads})
##' @rdname pars
##' @export
mkParInds <- function(parList) {
    if(!is.recursive(parList)) stop("parList must be a list")
    if(length(parList) == 1L) return(lapply(parList, seq_along))
    parInds <- mapply(`+`,
                      lapply(parList, seq_along),
                      c(0, cumsum(lapply(parList, length))[-length(parList)]),
                      SIMPLIFY = FALSE)
    names(parInds) <- names(parList) ## too paranoid?
    keepers <- sapply(parInds, length) > 0
    parInds[keepers]
}


initializeResp <- function(y, etastart = NULL, mustart = NULL,
                           offset = NULL, weights = NULL,
                           family){
    # taken mostly from mkRespMod
    
    ## y <- model.response(fr)
### Why consider there here?  They are handled in plsform.
    # offset <- model.offset(fr)
    # weights <- model.weights(fr)
    n <- length(y)
    etastart_update <- etastart
    if(length(dim(y)) == 1) {
    ## avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
### I really wish that the glm families in R were cleaned up.  All of
### this is such an ugly hack, just to handle one special case for the binomial
    rho <- new.env()
    rho$y <- if (is.null(y)) numeric(0) else y
    rho$etastart <- etastart
    rho$mustart <- mustart
    if (!is.null(offset)) {
        if (length(offset) == 1L) offset <- rep.int(offset, n)
        stopifnot(length(offset) == n)
        rho$offset <- unname(offset)
    } else rho$offset <- rep.int(0, n)
    if (!is.null(weights)) {
        stopifnot(length(weights) == n, all(weights >= 0))
        rho$weights <- unname(weights)
    } else rho$weights <- rep.int(1, n)
    stopifnot(inherits(family, "family"))
    rho$nobs <- n
    eval(family$initialize, rho)
    family$initialize <- NULL     # remove clutter from str output
    rho
}

fixFamily <- function(family) {
    if(!inherits(family, "family")) {
        if(inherits(family,"character")) {
            family <- get(family, mode = "function", envir = parent.frame())
        }
        if(inherits(family, "function")) {
            family <- family()
        }
    }
    return(family)
}

