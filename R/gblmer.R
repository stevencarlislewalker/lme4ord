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
##' @export
gblmer <- function(linFormula, bilinFormula,
                   data, family,
                   latentDims = 0L, loadingPen = 0L,
                   verbose = 0L, ...) {
    if(!any(inherits(data, "data.list"))) stop("data must be a data list")
    if(length(dim(data)) != 2L) stop("data list must be two dimensional")
    dIds <- names(dd <- dim(data))
    df <- as.data.frame(dims_to_vars(data))
    if(latentDims == 0L) {
        warning("no latent variables, calling glmer")
        return(glmer(linFormula, df, family, ...))
    }
    if(latentDims != 1L) stop("code for more than one latent variable not writen")
    # form1 <- formula
    # form2 <- as.formula(paste(". ~ 0 + (0 + latent |", names(dd[2], ")"))
    bilinFormula[[2]] <- linFormula[[2]] # use same response variables
    parsedLinFormula <- parsedForm <- glFormula(linFormula, df, family, ...)
    parsedBilinFormula <-               glFormula(bilinFormula, df, family, ...)
    reTrms <- joinReTrms(parsedBilinFormula$reTrms, parsedLinFormula$reTrms)
    parsedForm$reTrms <- reTrms

    theta <- reTrms$theta
    lower <- reTrms$lower
    dfun <- do.call(mkGlmerDevfun, parsedForm)
    rho <- environment(dfun)

                                        # which Zt@x elements
                                        # represent loadings
    rho$Zwhich <- rho$pp$Zt@i %in% (seq_len(dd[2]) - 1)
                                        # mapping from loadings to the
                                        # Zt@x elements that represent
                                        # loadings
    rho$Zind <- with(rho, pp$Zt@x[Zwhich])
    rho$Ztx <- rho$pp$Zt@x
    rho$loadInd <- 1:dd[1]

    dfunPrefix <- function(pars) {
        theta <- c(1, pars[-loadInd]) # the '1' is is for scalar
                                      # bilinear random effects
                                      # (FIXME: be more general)
        Ztx[Zwhich] <- pars[loadInd][Zind]
        trash <- pp$setZt(Ztx)
    }

    dfunSuffix <- function(pars) {
        p + loadingPen * sum(pars[loadInd]^2)
    }

    body(dfun) <- cBody(body(dfunPrefix), body(dfun), body(dfunSuffix))
    formals(dfun) <- setNames(formals(dfun), "pars")

    initLoadings <- svd(scale(t(data$Y)))$v[,1]
    
    #opt <- optim(c(initLoadings, theta[-1]), dfun, method = "L-BFGS-B",
    #             lower = c(rep(-Inf, dd[1]), lower),
    #             control = list(trace = 3))

                                        # here is the '-1' again for
                                        # scale bilinear random
                                        # effects
    opt <- lme4:::optwrap("bobyqa", dfun, c(initLoadings, theta[-1]), 
                   lower = c(rep(-Inf, dd[1]), lower), verbose = verbose)

    optLoadings <- opt$par[rho$loadInd]
    optTheta <- c(1, opt$par[-rho$loadInd])

    ## rho$control <- attr(opt,"control")
    rho$nAGQ <- 0
    opt$par <- optTheta

    body(dfun) <- body(dfun)[c(1, 5:9)]
    formals(dfun) <- setNames(formals(dfun), "theta")

    ## opt <- optimizeGlmer(dfun) # optimize without updating loadings
    ## opt$par <- optTheta
    dfun(optTheta)

    mer <- mkMerMod(environment(dfun), opt, parsedForm$reTrms, parsedForm$fr)

                                        # (FIXME: write specific
                                        # mkGlmerLatentMod)
    merList <- list()
    for(sl in names(getSlots("glmerMod"))) {
        merList[[sl]] <- slot(mer, sl)
    }
    do.call(new, c(list(Class = "gblmerMod"),
                   merList,
                   list(loadings = optLoadings)), quote = TRUE)
}

##' Class "gblmerMod"
##'
##' Class \code{"gblmerMod"} is san S4 class that extends
##' \code{"glmerMod"}
##' @name gblmerMod-class
##' @aliases gblmerMod-class
##' @docType class
##' @section Slots: in addition to the slots provided by
##' \code{"glmerMod"}, there is an additional \code{"loadings"} slot
##' containing the factor loadings
##' @keywords classes
##' @export
setClass("gblmerMod",
         representation(loadings = "numeric"),
         contains="glmerMod")

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
        # ------------------------------------------------------------
        # copied from mkReTrms (FIXME: break this out for reuse)
        flist <- as.data.frame(cbind(flist1, flist2))
        fnms <- names(flist)
        if (length(fnms) > length(ufn <- unique(fnms))) {
            flist <- flist[match(ufn, fnms)]
            asgn <- match(fnms, ufn)
        } else asgn <- seq_along(flist)
        names(flist) <- ufn
        flist <- do.call(data.frame, c(flist, check.names = FALSE))
        attr(flist, "assign") <- asgn
        # ------------------------------------------------------------
        list(Zt = rBind(Zt1, Zt2),
             theta = c(theta1, theta2),
             Lind = c(Lind1, Lind2 + max(Lind1)),
             Gp = c(Gp1, Gp2[-1] + max(Gp1)),
             lower = c(lower1, lower2),
             Lambdat = .bdiag(list(Lambdat1, Lambdat2)),
             flist = flist,
             cnms = c(cnms1, cnms2),
             Ztlist = c(Ztlist1, Ztlist2))
    })
}
