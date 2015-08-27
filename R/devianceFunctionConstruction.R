## ----------------------------------------------------------------------
## Deviance function construction
##
## used by:  structuredGLMM.R
## ----------------------------------------------------------------------


##' Make a general (Laplace) glmer deviance function
##'
##' @param y response vector
##' @param X fixed effects model matrix
##' @param Zt transposed random effects model matrix
##' @param Lambdat relative covariance factor
##' @param weights prior weights
##' @param offset known additive offsets
##' @param etastart,mustart starting values for linear predictor
##' (\code{etastart}) or mean (\code{mustart})
##' @param devfunEnv environment of the returned deviance function
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
##' @param mapToWeights function taking the \code{weigh} parameters
##' @param penFixef,penCovar,penLoads optional functions of each
##' parameter type, returning a scalar penalty for the deviance
##' function
##' @param family \code{\link{family}}
##' @param tolPwrss tolerance for penalized weighted residual sum of
##' squares
##' @param verbose verbose
##' @param maxit maximum number of PIRLS iterations
##' @param Lind optional indices mapping covariance parameters to the
##' non-zero elements of the relative covariance factor
##' @return a deviance function with an environment
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
                                 penFixef = NULL,
                                 penCovar = NULL,
                                 penLoads = NULL,
                                 family,
                                 tolPwrss = 1e-6,
                                 maxit = 30,
                                 verbose = 0L,
                                 Lind = NULL) {

    if(is.matrix(y)) stop("Only vector-valued responses allowed.\n",
                          "If you are trying to specify a binomial model,\n",
                          "please use weights to specify the total number\n",
                          "of Bernoulli trials.")

    ## silence no visible binding for global variable notes
    resp <- baseOffset <- lp0 <- setCovar <- pp <- setLoads <- setWeigh <-
        GQmat <- compDev <- fac <- NULL

    ## FIXME: institute lme4pureR possibility

    if(isLind <- !is.null(Lind)) {
        theta <- environment(mapToCovFact)$trans(initPars[parInds$covar])
    } else {
        theta <- Lambdat@x
        Lind <- seq_along(Lambdat@x)
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
    if(is.null(etastart))  etastart  <- family$linkfun(initializeEnv$mustart)

    devfunList <- list(Lind = Lind,
                       pp = lme4::merPredD$new(
                           X = X, Zt = Zt,
                           Lambdat = Lambdat,
                           Lind = Lind,
                           theta = as.double(theta),
                           n = nrow(X)),
                       resp = lme4::glmResp$new(
                           y = y, family = family,
                           weights = weights),
                       lp0 = etastart,
                       baseOffset = offset,
                       tolPwrss = tolPwrss,
                       maxit = maxit,
                       GQmat = GHrule(1), ## always Laplace approx
                       compDev = TRUE,
                       fac = NULL,
                       verbose = verbose,
                       setCovar = !is.null(parInds$covar),
                       setLoads = !is.null(parInds$loads),
                       setWeigh = !is.null(parInds$weigh),
                       setFixef = !is.null(parInds$fixef),
                       isLind = isLind,
                       mapToCovFact = mapToCovFact,
                       mapToModMat  = mapToModMat,
                       mapToWeights = mapToWeights,
                       parInds = parInds)


    devfun <- function(pars) {
        resp$setOffset(baseOffset)
        resp$updateMu(lp0)
        if(setCovar) {
            if(isLind) {
                pp$setTheta(as.double(environment(mapToCovFact)$trans(pars[parInds$covar])))
            } else {
                pp$setTheta(as.double(mapToCovFact(pars[parInds$covar])))
            }
        }
        if(setLoads) pp$setZt(as.double(mapToModMat(pars[parInds$loads])))
        if(setWeigh) resp$setWeights(as.double(mapToWeights(pars[parInds$weigh])))
        spars <- as.numeric(pars[parInds$fixef])
        offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
        resp$setOffset(offset)
        ## pp, resp, nAGQ, tol, maxit, verbose
        p <- lme4::glmerLaplaceHandle(pp$ptr(), resp$ptr(), 1, tolPwrss, maxit, verbose)
        ##p <- lme4:::glmerPwrssUpdate(pp, resp, tolPwrss, GQmat,
        ##                             compDev, fac, maxit, verbose)
        resp$updateWts()
        p
    }

    penFun <- mkPenaltyFun(penCovar, penFixef, penLoads, devfunEnv)
    if(!is.null(penFun)) {
        devfunList$penFun <- penFun
                                        # 12 means line 12 of the
                                        # devfun
        body(devfun)[[12]] <- quote(p + penFun(pars, parInds)) 
    }
    environment(devfun) <- list2env(devfunList, envir = devfunEnv)

                                        # initialize weights etc ...
    devfunEnv$resp$updateMu(etastart)
    devfunEnv$resp$updateWts()
    devfun(initPars)
    environment(devfun)$lp0 <- environment(devfun)$pp$linPred(1)

    return(devfun)
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

##' Make penalty function
##'
##' @param penCovar,penFixef,penLoads Penalty functions for
##' covariance, fixed effects, and loadings parameters.  Can be
##' \code{NULL}.
##' @param env Environment in which to evaluate the resulting
##' function.
##' @return A function with two arguments: 
##' \item{pars}{Parameter vector}
##' \item{parInds}{Named list of indices pointing to the elements of
##' \code{pars} corresponding to one of the three parameter types. The
##' names of this list must include \code{covar}, \code{fixef}, and/or
##' \code{loads} if the correponding type has a supplied
##' (i.e. non-null) penalty function.}
##' @export
##' @examples
##' penCovar <- function(x) x^2
##' penFixef <- NULL
##' penLoads <- function(gg) abs(gg) + 10
##' mkPenaltyFun(penCovar, penFixef, penLoads)(1:10, list(covar = 1:5, loads = 6:10))
mkPenaltyFun <- function(penCovar, penFixef, penLoads, env = environment()) {
    penFuncs <- list(penCovar, penFixef, penLoads)
    keepers <- !sapply(penFuncs, is.null)
    if(!any(keepers)) return(NULL)
    parTypes <- c("Covar", "Fixef", "Loads")
    parTypesName <- paste("pen", parTypes, sep = "")
    indsChar <- lapply(lapply(tolower(parTypes), c, list("$", "parInds")),
                       "[", c(2, 3, 1))
    parTypeCall <- lapply(lapply(rapply(indsChar, as.name, how = "list"), as.call),
                          function(inds) {
                              as.call(list(as.name("["), as.name("pars"), inds))})
    penSumCall <- c(list(quote(sum)),
                    mapply(call, parTypesName[keepers], parTypeCall[keepers]))
    env <- list2env(setNames(penFuncs, parTypesName), env)
    eval(call("function", as.pairlist(alist(pars = , parInds = )),
              as.call(unname(penSumCall))), env)
}

##' @param p exponent of the L-p norm
##' @param lambda multiplier for the penalty term
##' @rdname mkPenaltyFun
##' @export
mkPenLpNorm <- function(p = 2, lambda = 1) {
    if(p < 1) stop("p must be greater than or equal to 1")
    if(!(lambda > 0)) stop("lambda must be greater than 0")
    local({
        p <- p
        lambda <- lambda
        function(pars) lambda * sum(abs(pars^p))
    })
}

##' @param alpha,beta parameters of the generalized double Pareto
##' distribution (defaults from Murray et al 2013, JASA)
##' @rdname mkPenaltyFun
##' @export
mkPenPareto <- function(alpha = 3, beta = 1) {
    local({
        alpha <- alpha
        beta <- beta
        function(pars) 2 * (alpha + 1) * sum(log(1 + abs(pars)/beta))
    })
}

