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
##' \code{\link{lme4pureR}}?
##' @return a deviance function with an environment
##' @export
mkGeneralGlmerDevfun <- function(y, X, Zt, Lambdat,
                                 weights, offset,
                                 initPars, parInds,
                                 mapToCovFact, mapToModMat,
                                 family = binomial(),
                                 tolPwrss = 1e-6,
                                 verbose = 0L, pureR = FALSE) {
    if(pureR) {
        return(pirls(X, y, Zt, Lambdat, mapToCovFact,
                     initPars, weights, offset,
                     family = family, tol = tolPwrss))
    }
    devfun <- local({
        Lind <- seq_along(Lambdat@x)
        pp <- merPredD$new(X = X, Zt = Zt, Lambdat = Lambdat, Lind = Lind,
                           theta = as.double(Lambdat@x), n = nrow(X))
        resp <- glmResp$new(y = y, family = family, weights = weights)
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

##' @export
covar <- function(object, ...) UseMethod("covar")

##' @export
loads <- function(object, ...) loadings(object)

##' @export
pars <- function(object, ...) UseMethod("pars")

.covar <- function(pars, ind) pars[ind$covar]
.fixef <- function(pars, ind) pars[ind$fixef]
.loads <- function(pars, ind) pars[ind$loads]
