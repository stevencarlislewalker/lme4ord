##' Generalized bilinear mixed model
##'
##' @param linFormula a mixed model formula for the linear part of the
##' model
##' @param bilinFormula a mixed model formula for the bilinear part of
##' the model
##' @param data a \code{\link{data.list}} object
##' @param family a \code{\link{family}} object
##' @param latentDims number of latent dimensions in the bilinear
##' component of the model
##' @param loadingPen penalty for the size of the loadings
##' @param verbose passed to optimizer
##' @param ... arguments to be passed to \code{\link{glFormula}}
##' @importClassesFrom lme4 merMod glmerMod
##' @import multitable
##' @export
gblmer <- function(formula, data, family,
                   latentName = "latent", latentDims = 0L,
                   loadingsDim = 1L, loadingPen = 0L,
                   verbose = 0L, ...) {
    if(!any(inherits(data, "data.list"))) stop("data must be a data list")
    if(length(dim(data)) != 2L) stop("data list must be two dimensional")
    dIds <- names(dd <- dim(data))
    df <- as.data.frame(dims_to_vars(data))
    initGlmer <- glmer(formula, df, family, ...)
    if(latentDims == 0L) {
        warning("no latent variables, returning glmer results")
        return(initGlmer)
    }
    ## if(latentDims != 1L) stop("code for more than one latent variable not writen")
    if(loadingsDim == 2L) {
        datY <- data$Y
    } else if(loadingsDim == 1L) {
        datY <- t(data$Y)
    } else {
        stop("loadingsDim neight 1 nor 2")
    }
    U <- try(matrix(0, nrow = dd[loadingsDim], ncol = latentDims))
    if(inherits(U, "try-error")) stop("loadingsDim does not index a dimension of data")
    nFreeLoadings <- (dd[loadingsDim] * latentDims) - choose(latentDims, 2)
    U[lower.tri(U, TRUE)] <- 1:nFreeLoadings
    latentVarNames <- paste(latentName, 1:latentDims, sep = "")
    U <- setNames(as.data.frame(U), latentVarNames)
    latentData <- data.list(U, drop = FALSE, dimids = dIds[loadingsDim])
    data <- data + latentData
    df <- as.data.frame(dims_to_vars(data))

    # form1 <- formula
    # form2 <- as.formula(paste(". ~ 0 + (0 + latent |", names(dd[2], ")"))
    linFormula <- formula
    latentForm <- paste(latentVarNames, collapse = " + ")
    bilinFormula <- as.formula(paste(". ~ 0 + (0 + ", latentForm, " ||",
                                     dIds[-loadingsDim], ")"))
                                        # FIXME: replace -tv with more
                                        # general removal strategy
                                        # (e.g. with characters etc)
    bilinFormula[[2]] <- linFormula[[2]] # use same response variables

      parsedLinFormula <- parsedForm <- glFormula(  linFormula, df, family, ...)
    parsedBilinFormula <-               glFormula(bilinFormula, df, family, ...)
    reTrms <- joinReTrms(parsedBilinFormula$reTrms, parsedLinFormula$reTrms)
    parsedForm$reTrms <- reTrms

    theta <- reTrms$theta
    lower <- reTrms$lower
    dfun <- do.call(mkGlmerDevfun, parsedForm)
    dfun <- updateGlmerDevfun(dfun, reTrms)
    rho <- environment(dfun)

                                        # which Zt@x elements
                                        # represent loadings
    rho$Zwhich <- rho$pp$Zt@i %in% (seq_len(latentDims * dd[-loadingsDim]) - 1)
                                        # mapping from loadings to the
                                        # Zt@x elements that represent
                                        # loadings
    rho$Zind <- with(rho, pp$Zt@x[Zwhich])
    rho$Ztx <- rho$pp$Zt@x
    rho$loadInd <- 1:nFreeLoadings

    dfunPrefix <- function(pars) {
        loadings <- pars[loadInd]
        Ztx[Zwhich] <- loadings[Zind]
        trash <- pp$setZt(Ztx)
        pars <- c(rep(1, latentDims), pars[-loadInd])
                                        # the '1' is is for scalar
                                        # bilinear random effects
                                        # (FIXME: be more general)
    }

    dfunSuffix <- function(pars) {
        p + loadingPen * sum(loadings^2)
    }

    body(dfun) <- cBody(body(dfunPrefix), body(dfun), body(dfunSuffix))
    formals(dfun) <- setNames(formals(dfun), "pars")

    initLoadings <- svd(scale(datY))$v[, 1:latentDims, drop = FALSE]
    initLoadings <- initLoadings[lower.tri(initLoadings, TRUE)]
    
    #opt <- optim(c(initLoadings, theta[-1]), dfun, method = "L-BFGS-B",
    #             lower = c(rep(-Inf, dd[1]), lower),
    #             control = list(trace = 3))

                                        # here is the '1' again (with
                                        # a minus in front) for scale
                                        # bilinear random effects

    # initial parameters
    # order: (1) loadings, (2) covariance pars, (3) fixed effect pars
    parPointer <- setNames(cumsum(c(0, length(initLoadings),
                                    length(theta) - latentDims)) + 1,
                           c("loadings", "theta", "fixef"))
    initPars <- c(initLoadings, theta[-(1:latentDims)], rho$pp$beta(1))
    optLower <- c(rep(-Inf, length(initLoadings)), rho$lower[-(1:latentDims)])
                                        # optimize
    opt <- lme4:::optwrap("bobyqa", dfun, initPars, 
                          lower = optLower, verbose = verbose)

    optPar <- opt$par
    optNoLoadings <- c(rep(1, latentDims), opt$par[-rho$loadInd])
    optLoadings <- matrix(0, dd[loadingsDim], latentDims)
    optLoadings[lower.tri(optLoadings, TRUE)] <- optPar[rho$loadInd]
    dim(optLoadings) <- c(dd[loadingsDim], latentDims)
    colnames(optLoadings) <- latentVarNames
    rownames(optLoadings) <- dimnames(data)[[loadingsDim]]
    opt$par <- optNoLoadings    

    ## rho$control <- attr(opt,"control")
    # rho$nAGQ <- 0

    # body(dfun) <- body(dfun)[c(1, 5:9)]
    # formals(dfun) <- setNames(formals(dfun), "theta")

    ## opt <- optimizeGlmer(dfun) # optimize without updating loadings
    ## opt$par <- optTheta
    # dfun(optTheta)

    mer <- mkMerMod(environment(dfun), opt, parsedForm$reTrms, parsedForm$fr)

                                        # (FIXME: write specific
                                        # mkGlmerLatentMod)
    merList <- list()
    for(sl in names(getSlots("glmerMod"))) {
        merList[[sl]] <- slot(mer, sl)
    }
    do.call(new, c(list(Class = "gblmerMod"),
                   merList,
                   list(  loadings = optLoadings,
                            optPar = optPar,
                        parPointer = parPointer)), quote = TRUE)
}

##' Class "gblmerMod"
##'
##' Class \code{"gblmerMod"} is san S4 class that extends
##' \code{"glmerMod"}
##' @name gblmerMod-class
##' @aliases gblmerMod
##' @aliases gblmerMod-class
##' \section{Slots}{
##'   \describe{
##'     \item{\code{loadings}:}{factor loadings.}
##'   }
##' }
##' @keywords classes
##' @export
setClass("gblmerMod",
         representation(  loadings = "matrix",
                            optPar = "numeric",
                        parPointer = "numeric"),
         contains="glmerMod")

##' Variance-covariance matrix of the (co)variance parameters and
##' loadings (and sometimes of fixed effects)
##'
##' @param object a \code{\link{gblmerMod}} object
##' @param correlation should the correlation matrix be computed too?
##' @param ... not used
##' @export
vcov.gblmerMod <- function(object, correlation = TRUE, ...) {
                                        # calc.vcov.hess is a function
                                        # from vcov.merMod (FIXME:
                                        # breakout of vcov.merMod and
                                        # expose?)
    calc.vcov.hess <- function(h) {
	## ~= forceSymmetric(solve(h/2)[i,i]) : solve(h/2) = 2*solve(h)
        h <- tryCatch(solve(h),
                      error=function(e) matrix(NA,nrow=nrow(h),ncol=ncol(h)))
        ## i <- -seq_len(ntheta)
	## h <- h[i,i]
	forceSymmetric(h + t(h))
    }
    fixefIndices <- object@parPointer[3]:length(object@optPar)
    h <- object@optinfo$derivs$Hessian[fixefIndices, fixefIndices, drop = FALSE]
    V.hess <- calc.vcov.hess(h)
    
    bad.V.hess <- any(is.na(V.hess))
    if (!bad.V.hess) {
        e.hess <- eigen(V.hess,symmetric = TRUE,only.values = TRUE)$values
        if (min(e.hess) <= 0) bad.V.hess <- TRUE
    }
    if (!bad.V.hess) {
        V <- V.hess
    } else {
        stop("variance-covariance matrix computed ",
             "from finite-difference Hessian is\n",
             "not positive definite or contains NA values")
    }
    rr <- tryCatch(as(V, "dpoMatrix"), error = function(e)e)
    if (inherits(rr, "error")) {
	warning(gettextf("Computed variance-covariance matrix problem: %s;\nreturning NA matrix",
                         rr$message), domain = NA)
        rr <- matrix(NA,nrow(V),ncol(V))
    }
    if(correlation)
	rr@factors$correlation <-
	    as(rr, "corMatrix") else rr # (is NA anyway)
    rr
}


##' Refactor \code{loadings} into a generic function
##'
##' @param x object with loadings to extract
##' @param ... additional arguments
##' @export
loadings <- function(x, ...) {
    UseMethod("loadings")
}

##' @export
loadings.default <- function(x, ...) x$loadings

##' @export
loadings.gblmerMod <- function(x, ...) x@loadings

##' Concatenate the bodies of functions
##'
##' @param ... function bodies to combine
##' @export
cBody <- function(...) {
    l... <- list(...)
    l...[-1] <- lapply(l...[-1], "[", -1)
    l...asList <- lapply(l..., as.list)
    l...cList <- do.call(c, l...asList, quote = TRUE)
    as.call(l...cList)
}

##' Join two lists describing two sets of random effects terms
##'
##' @param reTrms1,reTrms2 results of \code{\link{glFormula}(.)$reTrms}
##' @export
joinReTrms <- function(reTrms1, reTrms2) {
    names(reTrms1) <- paste(names(reTrms1), 1, sep = "")
    names(reTrms2) <- paste(names(reTrms2), 2, sep = "")
    reTrms <- c(reTrms1, reTrms2)
    with(reTrms, {
        flist <- joinFlist(flist1, flist2)
        attr(flist, "assign") <- c(match(names(cnms1), names(flist)),
                                   match(names(cnms2), names(flist)))
        q <- c(nrow(Zt1), nrow(Zt2))
        nth <- c(length(theta1), length(theta2))
        nCnms <- c(length(cnms1), length(cnms2))
        list(Zt = rBind(Zt1, Zt2),
             theta = c(theta1, theta2),
             Lind = c(Lind1, Lind2 + max(Lind1)),
             Gp = c(Gp1, Gp2[-1] + max(Gp1)),
             lower = c(lower1, lower2),
             Lambdat = .bdiag(list(Lambdat1, Lambdat2)),
             flist = flist,
             cnms = c(cnms1, cnms2),
             Ztlist = c(Ztlist1, Ztlist2),
             q = q, nth = nth, nCnms = nCnms)
    })
}

## ------------------------------------------------------------
## adpated from mkReTrms
## ------------------------------------------------------------
joinFlist <- function(flist1, flist2) {
    flist <- c(flist1, flist2)
    fnms <- names(flist)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        flist <- flist[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    } else asgn <- seq_along(flist)
    names(flist) <- ufn
    flist <- do.call(data.frame, c(flist, check.names = FALSE))
    return(flist)
}

##' Get random effects terms from a fitted \code{merMod} object
##'
##' @param object \code{\link{merMod}} object
##' @return see \code{\link{mkReTrms}}
##' @export
getReTrms <- function(object, ...) {
    rts <- c("Zt", "Lambdat", "Lind", "theta", "lower", "flist", "cnms")
    getME(object, rts)
}
