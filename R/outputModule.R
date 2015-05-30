##' Structured GLMM class
##'
##' An \code{S3} class for generalized linear mixed models with
##' structured (co)variance terms, which are lists with the following
##' elements: \describe{ \item{opt}{Result of the optimizer (currently
##' \code{bobyqa} in the \code{minqa} package).}
##' \item{parsedForm}{Results of \code{\link{strucParseFormula}}.}
##' \item{dfun}{A function for computing the model deviance.  The
##' environment of \code{dfun} contains objects for representing the
##' model.}  \item{mc}{Matched call}}
##' @name strucGlmer-class
##' @rdname strucGlmer-class
##' @exportClass strucGlmer
setOldClass("strucGlmer")

##' @param opt,parsedForm,dfun,mc See desciption of \code{strucGlmer}
##' objects above
##' @rdname strucGlmer-class
##' @export
mkStrucGlmer <- function(opt, parsedForm, dfun, mc) {
    
    names(opt$par) <- names(parsedForm$initPars)

    ans <- list(opt = opt, parsedForm = parsedForm, dfun = dfun, mc = mc)
    class(ans) <- "strucGlmer"

                                        # update the initialized
                                        # parameters in the repeated
                                        # sparse matrices to optimized
                                        # values
    trash <- mapply(setInit, parsedForm$random,
                    covarPerTerm(ans),
                    loadsPerTerm(ans))
    ansRand <- try(lapply(ans$parsedForm$random, update), silent = TRUE)
    if(inherits(ansRand, "try-error")) {
        warning("couldn't reset initial values to estimated values,\n",
                "so model printout is probably misleading.\n",
                "perhaps check the object itself.")
    } else {
        ans$parsedForm$random <- ansRand
    }
    ans$parsedForm$Lambdat <- update(ans$parsedForm$Lambdat, covar(ans))
    ans$parsedForm$Zt      <- update(ans$parsedForm$Zt,      loads(ans))
    return(ans)
}

## ----------------------------------------------------------------------
## printing and summary
## ----------------------------------------------------------------------

strucGlmerMethTitle <- function() {
    paste("Structured GLMM",
                  "fit by maximum likelihood (Laplace Approx)",
                  collpase = " ")
}

##' @param x,object \code{strucGlmer} objects
##' @param digits number of significant digits
##' @param ... additional arguments to methods
##' @rdname strucGlmer-class
##' @method print strucGlmer
##' @export
print.strucGlmer <- function(x, digits = max(3, getOption("digits") - 3),  ...) {
    lme4:::.prt.methTit(strucGlmerMethTitle(), class(x))
    lme4:::.prt.family(lme4:::famlink(x, resp = x$parsedForm$devfunEnv$resp))
    lme4:::.prt.call(x$mc)
    
    llAIC <- getStrucLlikAIC(x)
    lme4:::.prt.aictab(llAIC$AICtab, 4)

    lapply(x$parsedForm$random, printReTrm)

    if(length(cf <- fixef(x)) > 0) {
	cat("Fixed Effects:\n")
	print.default(format(cf, digits = digits),
		      print.gap = 2L, quote = FALSE, ...)
    } else cat("No fixed effect coefficients\n")

    ## FIXME: optimizer warnings??
}

##' @param use.hessian use numerical hessian in covariance
##' calculations if available
##' @rdname strucGlmer-class
##' @method summary strucGlmer
##' @export
summary.strucGlmer <- function(object, use.hessian = TRUE, ...) {
    vc <- vcov(object, use.hessian = use.hessian)
    resp <- object$parsedForm$devfunEnv$resp
    famL <- lme4:::famlink(resp = resp)
    p <- length(coefs <- fixef(object))
    coefs <- cbind("Estimate" = coefs ,
                   "Std. Error" = sqrt(Matrix::diag(vc)))
    if (p > 0) {
        coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level = 0)
        colnames(coefs)[3] <- paste("z", "value")
        coefs <- cbind(coefs, "Pr(>|z)" =
                       2 * pnorm(abs(cf3), lower.tail = FALSE))
    }

    llAIC <- getStrucLlikAIC(object)

    structure(list(methTitle = strucGlmerMethTitle(),
                   objClass = class(object),
                   logLik = llAIC[["logLik"]],
                   family = famL$fami, link = famL$link,
                   coefficients = coefs,
                   residuals = residuals(object, type = "deviance"),
                   vcov = vc, varcor = VarCorr(object),
                   AICtab = llAIC[["AICtab"]], call = object$mc),
              class = "summary.strucGlmer")

}

##' @param correlation display correlation among fixed effects?
##' @param signif.stars display significance stars?
##' @param show.resids show residuals?
##' @rdname strucGlmer-class
##' @method print summary.strucGlmer
##' @export
print.summary.strucGlmer <- function(x, digits = max(3, getOption("digits") - 3),
                                     correlation = NULL, 
                                     signif.stars = getOption("show.signif.stars"),
                                     show.resids = TRUE, ...) {
    ## basically just stolen from lme4
    
    lme4:::.prt.methTit(x$methTitle, x$objClass)
    lme4:::.prt.family(x)
    lme4:::.prt.call(x$call); cat("\n")
    lme4:::.prt.aictab(x$AICtab); cat("\n")
    if(show.resids) {
        lme4:::.prt.resids(x$residuals, digits = digits)
    }
    lapply(x$parsedForm$random, printReTrm)
    p <- nrow(x$coefficients)
    if(p > 0) {
        cat("\nFixed effects:\n")
        printCoefmat(x$coefficients, zap.ind = 3,
                     digits = digits, signif.stars = signif.stars)
        if(is.null(correlation)) {
            correlation <- p <= 10
            if(!correlation) {
                nam <- deparse(substitute(x))
                if(length(nam) > 1 || nchar(nam) >= 32) nam <- "..."
                message(sprintf(paste(
                    "\nCorrelation matrix not shown by default, as p = %d > %d.",
                    "Use print(%s, correlation = TRUE) or",
                    "    vcov(%s)        if you need it\n", sep = "\n"),
                                p, 10, nam, nam))
            }
        } else if(!is.logical(correlation)) stop("'correlation' must be NULL or logical")
        if(correlation) {
            if(is.null(VC <- x$vcov)) VC <- vcov(x, correlation = TRUE)
            corF <- VC@factors$correlation
            if(is.null(corF)) {
                message("\nCorrelation of fixed effects could have been required in summary()")
                corF <- cov2cor(VC)
            }
            p <- ncol(corF)
            if(p > 1) {
                rn <- rownames(x$coefficients)
                rns <- abbreviate(rn, minlength = 11)
                cat("\nCorrelation of Fixed Effects:\n")
                corf <- matrix(format(round(corF@x, 3), nsmall = 3),
                               ncol = p,
                               dimnames = list(rns, abbreviate(rn, minlength = 6)))
                corf[!lower.tri(corf)] <- ""
                print(corf[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }               
}
    

## ----------------------------------------------------------------------
## fitted values, residuals, covariance, etc...
## ----------------------------------------------------------------------

##' @rdname strucGlmer-class
##' @export
residuals.strucGlmer <- function(object, ...) {
    r <- residuals(object$parsedForm$devfunEnv$resp, "deviance", ...)
    if (is.null(nm <- names(object$parsedForm$response))) nm <- seq_along(r)
    names(r) <- nm
    r
}

##' @param ranefTrms character vector naming random effects terms to
##' be used
##' @param fixef should fixed effects be used?
##' @rdname strucGlmer-class
##' @export
fitted.strucGlmer <- function(object, ranefTrms, fixef = TRUE, ...) {
    if(fixef) {
        fe <- as.numeric(object$parsedForm$fixed %*% fixef(object))
    } else fe <- 0

    if(missing(ranefTrms)) ranefTrms <- seq_along(object$parsedForm$random)
    if(!is.null(ranefTrms)) {
        re <- Reduce("+", ranef(object, type = "ZLu")[ranefTrms])
    } else re <- 0

    return(fe + re)
}

##' @importFrom Matrix forceSymmetric
calc.vcov.hess <- function(h) {
    ## ~= forceSymmetric(solve(h/2)[i,i]) : solve(h/2) = 2*solve(h)
    h <- tryCatch(solve(h),
                  error=function(e) matrix(NA,nrow=nrow(h),ncol=ncol(h)))
    forceSymmetric(h + t(h))
}

##' @param justFixef should matrices be over fixed effects only?
##' @rdname strucGlmer-class
##' @export
vcov.strucGlmer <- function(object, correlation = TRUE,
                            use.hessian = TRUE, justFixef = TRUE, ...) {
    if(use.hessian) {
        optPar <- pars(object)
        h <- try(lme4:::deriv12(object$dfun, optPar)$Hessian, silent = TRUE)
        hess.avail <- !inherits(h, "error-try")
        if(!hess.avail) {
            stop(shQuote("use.hessian"),
                 "=TRUE specified, ",
                 "but Hessian can't be computed")
        }
    } else {
        hess.avail <- FALSE
    }
    V <- object$parsedForm$devfunEnv$pp$unsc()
    if(hess.avail) {
        V.hess <- calc.vcov.hess(h)
        bad.V.hess <- any(is.na(V.hess))
        if(!bad.V.hess) {
            e.hess <- eigen(V.hess, symmetric = TRUE, only.values = TRUE)$values
            if(min(e.hess) <= 0) bad.V.hess <- TRUE
        }
    }
    if(use.hessian) {
        if(!bad.V.hess) {
            V <- V.hess
        } else {
            warning("variance-covariance matrix computed ",
                    "from finite-difference Hessian is\n",
                    "not positive definite or contains NA values: falling back to ",
                    "var-cov estimated from RX")
        }
    }

    if(justFixef && hess.avail && use.hessian) {
        i <- object$parsedForm$parInds$fixef
        V <- V[i, i, drop = FALSE]
    }
    
    rr <- tryCatch(as(V, "dpoMatrix"), error = function(e)e)
    if (inherits(rr, "error")) {
	warning(gettextf("Computed variance-covariance matrix problem: %s;\nreturning NA matrix",
                         rr$message), domain = NA)
        rr <- matrix(NA,nrow(V),ncol(V))
    }

    if(justFixef) {
        nmsX <- colnames(object$parsedForm$devfunEnv$pp$X)
    } else {
        nmsX <- names(pars(object))
    }
    dimnames(rr) <- list(nmsX, nmsX)

    if(correlation)
	rr@factors$correlation <- rr
	    ## if(!is.na(sigm)) as(rr, "corMatrix") else rr # (is NA anyway)
    rr
}

##' @param sigma,rdig for consistency
##' @importFrom nlme VarCorr
##' @rdname strucGlmer-class
##' @export
VarCorr.strucGlmer <- function(x, sigma = 1, rdig = 3) {
    lapply(x$parsedForm$random, VarCorr)
}

formatStrucVC <- function(varc) {
    whichNoNames <- sapply(lapply(varc, rownames), is.null)
    ans <- vector("list", length(varc))
    ans[whichNoNames] <- lme4:::formatVC(varc[whichNoNames])
}


getStrucLlikAIC <- function(object) {
    llik <- logLik(object)
    AICstats <- c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik), 
                  deviance = deviance(object), df.resid = df.residual(object))
    list(logLik = llik, AICtab = AICstats)
}

##' @importFrom lme4 isREML
##' @rdname strucGlmer-class
##' @export
isREML.strucGlmer <- function(x, ...) FALSE

##' @rdname strucGlmer-class
##' @export
df.residual.strucGlmer <- function(object, ...) {
    nobs(object) - length(pars(object))
}

##' @rdname strucGlmer-class
##' @export
nobs <- function(object, ...) nrow(object$parsedForm$fixed)

##' @importFrom stats deviance
##' @rdname strucGlmer-class
##' @export
deviance.strucGlmer <- function(object, ...) sum(residuals(object, ...)^2) ## object$opt$fval

##' @importFrom stats logLik
##' @rdname strucGlmer-class
##' @export
logLik.strucGlmer <- function(object, ...) {
    val = -0.5 * deviance(object, ...)
    nobs <- nobs(object)
    structure(val, nobs = nobs, nall = nobs, df = length(pars(object)),
              class = "logLik")
}

##' @importFrom stats formula
##' @rdname strucGlmer-class
##' @export
formula.strucGlmer <- function(x, ...) x$parsedForm$formula


## ----------------------------------------------------------------------
## parameter retrieval
## ----------------------------------------------------------------------

##' @rdname pars
##' @export
pars <- function(object, ...) UseMethod("pars")

##' Parameter retrieval for structured generalized linear mixed models
##'
##' @param object a \code{strucGlmer} fitted model object
##' @rdname pars
##' @export
pars.strucGlmer <- function(object, ...) object$opt$par

##' @rdname pars
##' @export
pars.glmerMod <- function(object, ...) unlist(getME(object, c("theta", "beta")))

.covar <- function(pars, ind) pars[ind$covar]
.fixef <- function(pars, ind) pars[ind$fixef]
.loads <- function(pars, ind) pars[ind$loads]

##' @param type character string giving the type of parameter
##' (e.g. \code{"fixef", "covar"})
##' @rdname pars
##' @export
getStrucGlmerPar <- function(object, type, ...) {
    parInds <- environment(object$dfun)$parInds
    optPar <- object$opt$par
    optPar[unlist(parInds[type])]
}

##' @rdname pars
##' @export
covar <- function(object, ...) UseMethod("covar")

##' @param ... not used
##' @rdname pars
##' @export
covar.strucGlmer <- function(object, ...) {
    unname(getStrucGlmerPar(object, "covar"))
}

##' @rdname pars
##' @export
loads <- function(object, ...) loadings(object)

##' @importFrom stats loadings
##' @rdname pars
##' @export
loadings.strucGlmer <- function(object, ...) {
    getStrucGlmerPar(object, "loads")
}

##' @importFrom nlme fixef 
##' @rdname pars
##' @export
fixef.strucGlmer <- function(object, ...) {
    setNames(getStrucGlmerPar(object, "fixef"),
             colnames(object$parsedForm$fixed))
}

##' @importFrom nlme ranef
##' @rdname pars
##' @export
ranef.strucGlmer <- function(object, type = c("u", "Lu", "ZLu"), ...) {
    type <- type[1]
    pp <- object$parsedForm$devfunEnv$pp
    structs <- object$parsedForm$random
    nms <- names(nRePerTrm <- environment(object$dfun)$nRePerTrm)
    if(type ==   "u") {
        re <- pp$u(1)
    } else {
        re <- pp$b(1)
        if(type == "ZLu") {
            b <- subRagByLens(re, nRePerTrm)
            multFn <- function(struc, re) as.numeric(as.matrix(t(struc$Zt), sparse = TRUE) %*% re)
            return(setNames(mapply(multFn, structs, b, SIMPLIFY = FALSE), nms))
        }
    }
    setNames(subRagByLens(re, nRePerTrm), nms)
}

##' @param nParPerTrm vector of the number of parameters per term
##' @param pars parameter vector (e.g. result of \code{covar} or
##' \code{loads}
##' @rdname pars
##' @export
parPerTerm <- function(nParPerTrm, pars) {
    if(is.null(pars)) pars <- numeric(0)
    whichThere <- (nParPerTrm > 0) & (!is.na(nParPerTrm))
    ans <- vector("list", length(whichThere))
    ans[whichThere] <- subRagByLens(pars, nParPerTrm[whichThere])
    names(ans) <- names(nParPerTrm)
    return(ans)
}

##' @rdname pars
##' @export
covarPerTerm <- function(object) {
    parPerTerm(environment(object$dfun)$nLambdatParPerTrm,
               covar(object))
}

##' @rdname pars
##' @export
loadsPerTerm <- function(object) {
    parPerTerm(environment(object$dfun)$nZtParPerTrm,
               loads(object))
}


## ##' Get components from a structured glmer model
## ##'
## ##' @param object a \code{\link{strucGlmer}} object
## ##' @param name a character vector naming the component
## ##' @export
## getSME <- function(object,
##                   name =
##                   c("y", "X", "Zt", "Lambdat",
##                     "")) {
## }

