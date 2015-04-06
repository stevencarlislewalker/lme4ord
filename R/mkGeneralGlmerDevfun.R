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
##' @return a deviance function with an environment
##' @importFrom lme4pureR pirls
##' @importFrom lme4 GHrule
##' @export
mkGeneralGlmerDevfun <- function(y, X, Zt, Lambdat,
                                 weights, offset,
                                 initPars, parInds,
                                 mapToCovFact, mapToModMat,
                                 family = binomial(),
                                 tolPwrss = 1e-6,
                                 verbose = 0L, pureR = FALSE,
                                 Lind = NULL) {
    if(pureR) {
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
        
    devfun <- local({
        Lind <- Lind
        pp <- lme4:::merPredD$new(X = X, Zt = Zt, Lambdat = Lambdat, Lind = Lind,
                           theta = as.double(theta), n = nrow(X))
        resp <- lme4:::glmResp$new(y = y, family = family, weights = weights)
        lp0 <- pp$linPred(1)
        baseOffset <- offset
        tolPwrss <- 1e-3
        GQmat <- GHrule(1)
        compDev <- TRUE
        fac <- NULL
        verbose <- 0L
        setCovar <- !is.null(parInds$covar)
        setLoads <- !is.null(parInds$loads)
        setFixef <- !is.null(parInds$fixef)
        function(pars) {
            resp$setOffset(baseOffset)
            resp$updateMu(lp0)
            if(setCovar) pp$setTheta(as.double(mapToCovFact(pars[parInds$covar])))
            if(setLoads) pp$setZt(as.double(mapToModMat(pars[parInds$loads])))
            spars <- as.numeric(pars[parInds$fixef])
            offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
            resp$setOffset(offset)
            p <- lme4:::glmerPwrssUpdate(pp, resp, tolPwrss, GQmat,
                                         compDev, fac, verbose)
            resp$updateWts()
            p
        }
    })
    rho <- environment(devfun)
    rho$mapToCovFact <- mapToCovFact
    rho$mapToModMat <- mapToModMat
    rho$parInds <- parInds

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
