## ----------------------------------------------------------------------
## Output module
##
## used by:  strucGlmer.R
## uses:  reTrmStructs.R, repSparse.R
## ----------------------------------------------------------------------


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
    
    if(optSuccess <- !inherits(opt, "try-error")) names(opt$par) <- names(parsedForm$initPars)

    ans <- list(opt = opt, parsedForm = parsedForm, dfun = dfun, mc = mc)
    class(ans) <- "strucGlmer"

    if(!optSuccess) return(ans)

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
    .prt.methTit(strucGlmerMethTitle(), class(x))
    .prt.family(famlink(x, resp = x$parsedForm$devfunEnv$resp))
    .prt.call(x$mc)

    if(inherits(x$opt, "try-error")) {
        cat("Model failed to converge with error:\n")
        cat(x$opt)
        return(invisible())
    }
    
    llAIC <- getStrucLlikAIC(x)
    .prt.aictab(llAIC$AICtab, 4)

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
    famL <- famlink(resp = resp)
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
    ## basically just stolen from lme4, but with tweaks
    
    .prt.methTit(x$methTitle, x$objClass)
    .prt.family(x)
    .prt.call(x$call); cat("\n")
    .prt.aictab(x$AICtab); cat("\n")
    if(show.resids) {
        .prt.resids(x$residuals, digits = digits)
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
##' be used (if \code{NULL} no terms are used)
##' @param fixef should fixed effects be used?
##' @rdname strucGlmer-class
##' @export
fitted.strucGlmer <- function(object, ranefTrms, fixef = TRUE, ...) {
    if(fixef) {
        fe <- as.numeric(model.matrix(object) %*% fixef(object))
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
        h <- try(deriv12(object$dfun, optPar)$Hessian, silent = TRUE)
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
    ans[whichNoNames] <- formatVC(varc[whichNoNames])
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

##' @importFrom stats deviance
##' @rdname strucGlmer-class
##' @export
deviance.strucGlmer <- function(object, ...) sum(residuals(object, ...)^2) ## object$opt$fval

##' @importFrom stats logLik
##' @rdname strucGlmer-class
##' @export
logLik.strucGlmer <- function(object, ...) {
    val = laplace(object, ...)
    nobs <- nobs(object)
    structure(val, nobs = nobs, nall = nobs, df = length(pars(object)),
              class = "logLik")
}

##' @importFrom stats formula
##' @rdname strucGlmer-class
##' @export
formula.strucGlmer <- function(x, ...) x$parsedForm$formula

##' @importFrom stats family
##' @rdname strucGlmer-class
##' @export
family.strucGlmer <- function(object, ...) object$parsedForm$devfunEnv$resp$family

##' @importFrom stats weights
##' @rdname strucGlmer-class
##' @export
weights.strucGlmer <- function(object, ...) object$parsedForm$devfunEnv$resp$weights

##' @rdname strucGlmer-class
##' @export
getOffset <- function(object) object$parsedForm$devfunEnv$baseOffset

##' @importFrom stats model.matrix
##' @rdname strucGlmer-class
##' @export
model.matrix.strucGlmer <- function(object, ...) model.matrix(object$parsedForm)

## ----------------------------------------------------------------------
## simulations
## ----------------------------------------------------------------------

##' @param nsim number of simulations
##' @param seed random seed (see \code{\link{set.seed}})
##' @importFrom stats simulate
##' @rdname strucGlmer-class
##' @export
simulate.strucGlmer <- function(object, nsim = 1, seed = NULL, ...) {
    replicate(nsim, {
        fe <- getOffset(object) + fitted(object, ranefTrms = NULL)
        re <- Reduce("+", lapply(getReTrm(object), simReTrm)) ## Reduce is safer than sapply
        fam <- family(object)
        familySimFun(fam)(weights(object),
                          nobs(object),
                          fam$linkinv(fe + re))
    })
}

##' Simulate from a parsed formula
##' 
##' @param nsim number of simulations (FIXME: inconsistent between
##' binomial and poisson)
##' @param parsedForm result of \code{strucParseFormula}
##' @param family family object
##' @param weights optional weights
##' @export
simStrucParsedForm <- function(parsedForm, family = binomial,
                               weights, nsim) {
                               ## ranefTrms, fixef = TRUE) {
    if(missing(weights)) weights <- rep(1, length(parsedForm$response))
    with(parsedForm, {
        reMM <- t(as(Zt, "dgCMatrix")) %*% t(as(Lambdat, "dgCMatrix"))
        feMM <- fixed
        fe <- as.numeric(feMM %*% initPars[parInds$fixef])
        re <- as.numeric(reMM %*% rnorm(ncol(reMM)))
        simFun <- simfunList[[family()$family]]
        return(simFun(weights, nsim, family()$linkinv(fe + re)))
    })
}

## ----------------------------------------------------------------------
## print random effects terms
## ----------------------------------------------------------------------

##' Print random effects term
##' 
##' @param object \code{\link{repSparse}} object
##' @param forSummary print for \code{\link{summary}} instead of
##' \code{\link{print}}?
##' @param ... additional arguments
##' @export
printReTrm <- function(object, forSummary = FALSE, ...) {
    UseMethod("printReTrm")
}

##' @param x \code{reTrmStruct} object
##' @rdname printReTrm
##' @aliases print.reTrmStruct
##' @method print reTrmStruct
##' @export
print.reTrmStruct <- function(x, ...) {
    objectsFromSetting <- c("Zt", "Lambdat",
                            "lowerLoads", "upperLoads",
                            "lowerCovar", "upperCovar")
    if(!all(objectsFromSetting %in% names(x))){
        .title <- paste("Random effects term (class: ", class(x)[1], "):", sep = "")
        cat(.title, "\n")
        cat("Structure has not fully been set --",
            "use setReTrm\n")
        return(invisible(x))
    }
    printReTrm(x, ...)
}

.printPars <- function(description = "parameters: ", value) {
    if((length(value) > 0L) && (!is.na(value)) && (!is.null(value))) {
        cat(description, value, "\n")
    }
}

.printVC <- function(description = "variance-correlation: ", value) {
    if((nrow(value[[1]]) < 6) && (!is.null(rownames(value[[1]])))) {
        cat(description, "\n")
        print(formatVC(value), quote = FALSE, digits = 3)
    }
}

##' @rdname printReTrm
##' @export
printReTrm.default <- function(object, forSummary = FALSE, ...) {
    .title <- paste("Random effects term (class: ", class(object)[1], "):", sep = "")
    cat (.title, "\n")
    cpd <- if(length(cp <- getInit(object$Lambdat)) > 1L) {
        "  covariance parameters: " } else "  covariance parameter:  "
    lpd <- if(length(lp <- getInit(object$Zt)) > 1L) {
        "  loadings parameters:   " } else "  loadings parameter:    "
    vc <- structure(list(VarCorr(object)), names = object$grpName, useSc = FALSE)
    vcd <- "  variance-correlation:  "
    .printPars(cpd, cp)
    .printPars(lpd, lp)
    .printVC  (vcd, vc)
}

##' @rdname printReTrm
##' @export
printReTrm.factAnal <- function(object, forSummary = FALSE, ...) {
    .title <- paste("Random effects term (class: ", class(object)[1], "):", sep = "")
    cat (.title, "\n")
    trans <- environment(object$Zt$trans)$Btrans
    trmDims <- environment(trans)$trmDims
    loadMat <- matrix(trans(getInit(object$Zt)),
                      trmDims["nVar"], trmDims["nAxes"])
    
    .trmDims <- format(trmDims)
    latentAxes <- if(trmDims["nAxes"] > 1) " latent axes\n" else " latent axis\n"
    cat ("",
         .trmDims["nObs"],  " multivariate observations\n",
         .trmDims["nVar"],  " variables\n",
         .trmDims["nAxes"], latentAxes)

    loadSums <- cbind(mean   = format(apply(loadMat, 2, mean)),
                      stdDev = format(apply(loadMat, 2, sd)))
    rownames(loadSums) <- paste("axis", format(1:trmDims["nAxes"]))
    print(as.table(loadSums))
}

##' @rdname printReTrm
##' @export
printReTrm.nlmeCorStruct <- function(object, forSummary = FALSE, ...) {
    .title <- paste("Random effects term (class: ", class(object)[1], "):", sep = "")
    cat (.title, "\n")
    cpd <- "  corStruct object: "
    transEnv <- environment(object$Lambdat$trans)
    corObj <- transEnv$object
    vc <- structure(list(VarCorr(object)), names = object$grpName, useSc = FALSE)
    vcd <- "  variance-correlation:  "
    print(corObj)
    if(transEnv$sigExists) {
        cat("  Standard deviation multiplier: ", getInit(object$Lambdat), "\n")
    }
    .printVC(vcd, vc)
}

## FIXME: write more specific printReTrm methods for different classes


## ----------------------------------------------------------------------
## retrieval
## ----------------------------------------------------------------------


##' Get random effects structures from a strucGlmer object
##'
##' @param object \code{\link{strucGlmer}} or
##' \code{\link{strucParseFormula}} object
##' @param name Names of the random effects structures (if missing, a
##' list of the names are returned).  The naming convention for
##' \code{reTrmStruct} objects is
##' \code{grpFacName.reTrmStructClass}. Partial matching of names is
##' allowed.
##' @param drop If only one term is selected, return the
##' \code{reTrmStruct} objectitself rather than a length-one list with
##' the object.
##' @return a fitted random effects structure
##' @seealso \code{\link{setReTrm}}
##' @export
getReTrm <- function(object, name, drop = TRUE) {
    if(inherits(object, "strucGlmer")) object <- object$parsedForm
    structs <- object$random
    if(missing(name)) {
        #message("available reTrmStruct objects:\n",
        #        paste(names(structs), collapse = ", "))
        #return(invisible())
        return(structs)
    }
    nameInd <- which(sapply(lapply(names(structs), "==", name), any))
    if(length(nameInd) == 0L) {
        stop("could not find matching struct. ",
             "please try one of the following:\n",
             paste(names(structs), collapse = ", "))
    }
    if(drop & (length(nameInd) == 1L)) {
        return(structs[[nameInd]])
    } else {
        return(structs [nameInd] )
    }
}


## ----------------------------------------------------------------------
## Compression
## ----------------------------------------------------------------------

##' Compress strucGlmer object
##'
##' This feature is not yet written, but is important enough that
##' current design decisions should reflect it's future possibility.
##' This proposal will help reduce the costs of keeping redundant
##' information in \code{\link{strucGlmer}} objects. The idea is to
##' store various components on disk, and then only retreive them when
##' necessary. This will require very consistent use of extractor
##' functions that have versions that can take a filename or
##' connection and retreive the component that way.
##'
##' @section What kind of redundancy are we talking about?:
##' \code{\link{strucGlmer}} objects contain a
##' \code{\link{strucParseFormula}} object and a deviance
##' function. The environment of this deviance function contains many
##' objects, some of which are essentially replicated in the parsed
##' formula object. In particular, the matrices in the parsed formula
##' are \code{\link{repSparse}} objects whereas the matrices in the
##' environment of the deviance function are \code{Matrix} package
##' objects linked to \code{C++} objects through external pointers.
##' 
##' @section Why do we want redundancy?:
##' This redundancy makes the output module code easier to maintain,
##' while retaining a fast optimization module. In partiular, the
##' \code{\link{repSparse}} objects in the parsed formula make it
##' easier to write consistent and resuseable code, whereas the
##' \code{Matrix}/\code{C++} objects in the environment of the
##' deviance function facilitate relatively faster linear
##' algebra.
##'
##' @param object a \code{\link{strucGlmer}} object.
##' @param components components to compress.
##' @param ... not used yet.
##' @return an error message pointing to this help page.
##' @export
compressStrucGlmer <- function(object, components, ...) {
    stop("Compression not yet written.\n",
         "Please see ?compressStrucGlmer for more info.")
}

