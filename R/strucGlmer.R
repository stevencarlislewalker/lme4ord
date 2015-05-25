##' Generalized linear mixed model with structured (co)variance terms
##'
##' Fits a generalized linear mixed model with structured (co)variance
##' terms by the following fully modularized steps:
##' \describe{
##' \item{\code{\link{strucParseFormula}}}{parses the mixed model
##' formula with structured terms}
##' \item{\code{\link{mkGeneralGlmerDevfun}}}{constructs a generalized
##' linear mixed model deviance function}
##' \item{\code{\link{bobyqa}}}{optimizer the deviance function}
##' \item{\code{\link{mkStrucGlmer}}}{constructs an object of
##' \code{\link{strucGlmer-class}}}}
##'
##' @param formula extended mixed model formula
##' @param data data frame
##' @param addArgs list of additional arguments to pass to
##' @param optVerb verbose
##' \code{\link{setReTrm}} methods
##' @param ... further arguments to \code{\link{mkGeneralGlmerDevfun}}
##' @export
strucGlmer <- function(formula, data, family, addArgs = list(), optVerb = 0L,
                       weights = NULL, offset = NULL, etastart = NULL,
                       devfunOnly = FALSE,
                         ...) {

    mc <- match.call()
    
    cat("\nConstructing vectors and matrices...\n")
    parsedForm <- strucParseFormula(formula, data, addArgs, ...)

    cat("\nConstructing deviance function...\n")
    dfun <- mkGeneralGlmerDevfun(y = parsedForm$response,
                                 X = parsedForm$fixed,
                                      Zt = as(parsedForm$Zt,      "dgCMatrix"),
                                 Lambdat = as(parsedForm$Lambdat, "dgCMatrix"),
                                 weights = weights,
                                 offset = offset,
                                 etastart = etastart,
                                 initPars = parsedForm$initPars,
                                 parInds = parsedForm$parInds,
                                 mapToCovFact = mkSparseTrans(parsedForm$Lambdat),
                                 mapToModMat = mkSparseTrans(parsedForm$Zt),
                                 devfunEnv = parsedForm$devfunEnv,
                                 family = family,
                                 ...)
    if(devfunOnly) return(dfun)

    cat("\nOptimizing deviance function...\n")
    opt <- minqa::bobyqa(parsedForm$initPars, dfun,
                         lower = parsedForm$lower,
                         upper = parsedForm$upper,
                         control =
                         list(iprint = optVerb,
                              rhobeg = 0.0002,
                              rhoend = 2e-7))
    
    cat("\nPreparing output...\n")
    mkStrucGlmer(opt, parsedForm, dfun, mc)
}

##' Structured GLMM class
##'
##' An \code{S3} class for generalized linear mixed models with
##' structured (co)variance terms, which are lists with the following
##' elements:
##' \describe{
##'   \item{opt}{Result of the optimizer (currently \code{\link{minqa}}).}
##'   \item{parsedForm}{Results of \code{\link{strucParseFormula}}.}
##'   \item{dfun}{A function for computing the model deviance.
##'               The environment of \code{dfun} contains objects for
##'               representing the model.}
##'   \item{mc}{Matched call}}
##' @name strucGlmer-class
##' @rdname strucGlmer-class
##' @exportClass strucGlmer
setOldClass("strucGlmer")

## ----------------------------------------------------------------------
## printing and summary
## ----------------------------------------------------------------------

strucGlmerMethTitle <- function() {
    paste("Structured GLMM",
                  "fit by maximum likelihood (Laplace Approx)",
                  collpase = " ")
}

##' @param x,object \code{strucGlmer} objects
##' @param ... additional arguments to methods
##' @rdname strucGlmer-class
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

##' @rdname strucGlmer-class
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

##' @rdname strucGlmer-class
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
    

##' @rdname strucGlmer
##' @export
residuals.strucGlmer <- function(object, ...) {
    r <- residuals(object$parsedForm$devfunEnv$resp, "deviance", ...)
    if (is.null(nm <- names(object$parsedForm$response))) nm <- seq_along(r)
    names(r) <- nm
    ## if (!is.null(na.action <- attr(model.frame(object), "na.action")))
    ##     r <- naresid(na.action, r)
    r
}

##' @rdname strucGlmer
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

##' @rdname strucGlmer
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

##' @importFrom nlme VarCorr
##' @rdname strucGlmer
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
deviance.strucGlmer <- function(object, ...) object$opt$fval

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


##' @param type character string giving the type of parameter
##' (e.g. \code{"fixef", "covar"})
##' @rdname strucGlmer
##' @export
getStrucGlmerPar <- function(object, type, ...) {
    parInds <- environment(object$dfun)$parInds
    optPar <- object$opt$par
    optPar[unlist(parInds[type])]
}


##' @importFrom nlme fixef 
##' @rdname strucGlmer
##' @export
fixef.strucGlmer <- function(object, ...) {
    setNames(getStrucGlmerPar(object, "fixef"),
             colnames(object$parsedForm$fixed))
}

##' @rdname strucGlmer
##' @export
covar.strucGlmer <- function(object, ...) {
    unname(getStrucGlmerPar(object, "covar"))
}

##' @importFrom stats loadings
##' @rdname strucGlmer
##' @export
loadings.strucGlmer <- function(object, ...) {
    getStrucGlmerPar(object, "loads")
}

##' @rdname strucGlmer
##' @export
pars.strucGlmer <- function(object, ...) object$opt$par

subRagByLens <- function(x, lens) {
    split(x, rep(seq_along(lens), lens)) ## split no good ... order of levels !
}

##' @importFrom nlme ranef
##' @rdname strucGlmer-class
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

strucGlmerParPerTerm <- function(lens, pars) {
    if(is.null(pars)) pars <- numeric(0)
    whichThere <- (lens > 0) & (!is.na(lens))
    ans <- vector("list", length(whichThere))
    ans[whichThere] <- subRagByLens(pars, lens[whichThere])
    names(ans) <- names(lens)
    return(ans)
}

##' @rdname strucGlmer
##' @export
covarPerTerm <- function(object) {
    strucGlmerParPerTerm(environment(object$dfun)$nLambdatParPerTrm,
                         covar(object))
}

##' @rdname strucGlmer
##' @export
loadsPerTerm <- function(object) {
    strucGlmerParPerTerm(environment(object$dfun)$nZtParPerTrm,
                         loads(object))
}

##' Construct random effects structures
##'
##' Construct random effects structures from a formula and data.
##'
##' @param splitFormula results of \code{\link{splitForm}} (or a
##' \code{formula} itself)
##' @param data data
##' @seealso \code{\link{setReTrm}} for initializing the structures as
##' is required for model fitting, and \code{\link{splitForm}} for
##' breaking the formula up into terms.
##' @return A list of random effects structures, each associated with
##' the terms in a generalized mixed model formula of roughly the
##' following form: \code{(linForm1 | grpFac1) + (linForm2 | grpFac2)
##' + ... + specialStruc(linForm3 | grpFac3) + ...}.  Each such
##' structure may be unstructured with class `unstruc` (indicating
##' standard `lme4` structures, e.g. the first two terms above) or
##' class of same name as the function in the formula
##' (e.g. \code{specialStruc}).  Each structure also contains (at a
##' minimum) the following elements (although \code{\link{setReTrm}}
##' adds further structure):
##' 
##' \item{modMat}{A raw model matrix for which the structure is
##' defined.  This matrix is obtained by evaluating
##' \code{model.matrix(linForm, data)}, where \code{linForm | grpFac}
##' is the random effects formula defining the term.}
##'
##' \item{grpFac}{If present, the grouping factor associated with the
##' structure, \code{NA} otherwise.}
##'
##' \item{grpName}{If present, the name of the grouping factor,
##' \code{NA} otherwise.}  \item{addArgs}{If present, any additional
##' arguments passed to the special function in the formula.}
##' @export
mkReTrmStructs <- function(splitFormula, data) {
    if(inherits(splitFormula, "formula")) splitFormula <- splitForm(splitFormula)
    reTrmsList <- lapply(splitFormula$reTrmFormulas,
                         getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- paste(sapply(reTrmsList, "[[", "grpName"),
                               splitFormula$reTrmClasses, sep = ".")
    nUnStr <- sum(splitFormula$reTrmClasses == "unstruc")
    for(i in seq_along(reTrmsList)) {
        clsi <- splitFormula$reTrmClasses[[i]]
        if(clsi != "unstruc") reTrmsList[[i]]$addArgs <- splitFormula$reTrmAddArgs[[i - nUnStr]]
        class(reTrmsList[[i]]) <- c(clsi, "reTrmStruct")
    }
    return(reTrmsList)
}

##' Parse a mixed model formula
##'
##' @param formula mixed model formula
##' @param data an object coercible to data frame
##' @param addArgs list of additional arguments to
##' \code{\link{setReTrm}} methods
##' @param reTrmsList if \code{NULL} \code{\link{mkReTrmStructs}} is used
##' @param ... additional parameters to \code{\link{as.data.frame}}
##' @return A list with components:
##' \item{response}{The response vector}
##' \item{fixed}{The fixed effects model matrix}
##' \item{random}{List of random effects terms returned by
##' \code{\link{mkReTrmStructs}} and updated with
##' \code{\link{setReTrm}}.}
##' \item{Zt}{The sparse transposed random effects model matrix, of
##' class \code{dgCMatrix}.}
##' \item{Lambdat}{The sparse transposed relative covariance factor,
##' of class \code{dgCMatrix}}
##' \item{initPars}{The initial parameter vector}
##' \item{parInds}{List of indices to the different types of
##' parameters.}
##' \item{lower, upper}{Vectors of lower and upper bounds on
##' \code{initPars}.}
##' \item{devfunEnv}{Environment of the deviance function.}
##' \item{formula}{Model formula.}
##' @rdname strucParseFormula
##' @export
strucParseFormula <- function(formula, data, addArgs = list(), reTrmsList = NULL, ...) {
                                        # get and construct basic
                                        # information: (1) list of
                                        # formulas, (2) data, and (3)
                                        # the environment that will
                                        # eventually become the
                                        # environment of the deviance
                                        # function (i try to keep a
                                        # reference to this
                                        # environment in lots of
                                        # places)
    sf        <- splitForm(formula)
    data      <- as.data.frame(data, ...)
    devfunEnv <- new.env()

                                        # initially set up the random
                                        # effects term structures
    if(is.null(reTrmsList)) reTrmsList <- mkReTrmStructs(sf, data)

                                        # extract the respose, fixed
                                        # effects model matrix, and
                                        # list of random effects
                                        # structures
    response <- model.response(model.frame(sf$fixedFormula, data))
    fixed    <- model.matrix(sf$fixedFormula, data)
    random   <- lapply(reTrmsList, setReTrm, addArgs = addArgs, devfunEnv = devfunEnv)

                                        # lists of repeated sparse
                                        # matrices
    ZtList      <- lapply(random, "[[",      "Zt")
    LambdatList <- lapply(random, "[[", "Lambdat")

                                        # bind the lists together
    ZtBind      <- .bind(     ZtList, "row")
    LambdatBind <- .bind(LambdatList, "diag")

                                        # ensure that the order is
                                        # appropriate for coercing to
                                        # dgCMatrix objects
    Zt      <- standardSort(     ZtBind)
    Lambdat <- standardSort(LambdatBind)

                                        # get initial values for model
                                        # parameters
    init <- list(covar = getInit(Lambdat),
                 fixef = rep(0, ncol(fixed)),
                 loads = getInit(Zt))
    parInds <- mkParInds(init)
    initPars <- unlist(init)

                                        # get lower and upper bounds
                                        # on parameters for the
                                        # optimizer
    lowerLoadsList <- lapply(random, "[[", "lowerLoads")
    upperLoadsList <- lapply(random, "[[", "upperLoads")
    lowerCovarList <- lapply(random, "[[", "lowerCovar")
    upperCovarList <- lapply(random, "[[", "upperCovar")
    lower <- c(unlist(lowerCovarList),
               rep(-Inf, ncol(fixed)),
               unlist(lowerLoadsList))
    upper <- c(unlist(upperLoadsList),
               rep( Inf, ncol(fixed)),
               unlist(upperCovarList))
    names(lower) <- names(upper) <- names(initPars)

                                        # fill the environment of the
                                        # deviance function with those
                                        # objects that depend on the
                                        # order of random effects
                                        # terms
    devfunEnv <- list2env(list(nRePerTrm = sapply(LambdatList, nrow),
                               nLambdatParPerTrm = sapply(LambdatList, parLength),
                                    nZtParPerTrm = sapply(     ZtList, parLength),
                               reTrmClasses = sf$reTrmClasses),
                          envir = devfunEnv)

                                        # fill the environments of the
                                        # transformation functions
                                        # with objects that depend on
                                        # the order of random effects
                                        # terms
    random <- lapply(random, update)

    return(list(response = response, fixed = fixed, random = random,
                Zt = Zt, Lambdat = Lambdat,
                initPars = initPars, parInds = parInds,
                lower = lower, upper = upper,
                devfunEnv = devfunEnv,
                formula = formula))
}


##' @param parsedForm result of \code{strucParseFormula}
##' @param family family object
##' @param weights optional weights
##' @rdname strucParseFormula
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


##' Split a formula
##'
##' @param formula Generalized mixed model formula
##' @rdname splitForm
##' @export
splitForm <- function(formula) {

    specials <- findReTrmClasses()
                                        # ignore any specials not in
                                        # formula
    specialsToKeep <- sapply(lapply(specials, grep,
                                    x = as.character(formula[[length(formula)]])), length) > 0L
    specials <- specials[specialsToKeep]

    ## Recursive function: (f)ind (b)ars (a)nd (s)pecials
    ## cf. fb function in findbars (i.e. this is a little DRY)
    fbas <- function(term) {
        if (is.name(term) || !is.language(term)) return(NULL)
        for (sp in specials) if (term[[1]] == as.name(sp)) return(term)
        if (term[[1]] == as.name("(")) return(term)
        stopifnot(is.call(term))
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(fbas(term[[2]]))
        c(fbas(term[[2]]), fbas(term[[3]]))
    }
    formula <- expandDoubleVerts(formula)
                                        # split formula into separate
                                        # random effects terms
                                        # (including special terms)
    formSplits <- fbas(formula)
                                        # vector to identify what
                                        # special (by name), or give
                                        # "(" for standard terms, or
                                        # give "|" for specials
                                        # without a setReTrm method
    formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
    as.character(formSplits[[1]])
                                        # warn about terms without a
                                        # setReTrm method
    badTrms <- formSplitID == "|"
    if(any(badTrms)) {
        stop("can't find setReTrm method(s)\n",
             "use findReTrmClasses() for available methods")
        # FIXME: coerce bad terms to unstructured as attempted below
        warning(paste("can't find setReTrm method(s) for term number(s)",
                      paste(which(badTrms), collapse = ", "),
                      "\ntreating those terms as unstructured"))
        formSplitID[badTrms] <- "("
        fixBadTrm <- function(formSplit) {
            as.formula(paste(c("~(", as.character(formSplit)[c(2, 1, 3)], ")"),
                             collapse = " "))[[2]]
        }
        formSplits[badTrms] <- lapply(formSplits[badTrms], fixBadTrm)
    }

                                        # capture additional arguments
    reTrmAddArgs <- lapply(formSplits, "[", -2)[!(formSplitID == "(")]
                                        # remove these additional
                                        # arguments
    formSplits <- lapply(formSplits, "[", 1:2)
                                        # standard RE terms
    formSplitStan <- formSplits[formSplitID == "("]
                                        # structured RE terms
    formSplitSpec <- formSplits[!(formSplitID == "(")]

    if(length(formSplitSpec) == 0) stop(
                 "no special covariance structures. ",
                 "please use lmer or glmer")


    fixedFormula <- formula(paste(formula[[2]], "~",
                                  as.character(noSpecials(nobars(formula)))[[3]]))
    reTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                       lapply(formSplitSpec, "[[", 2))
    reTrmClasses <- c(rep("unstruc", length(formSplitStan)),
                      sapply(lapply(formSplitSpec, "[[", 1), as.character))
    
    return(list(fixedFormula  = fixedFormula,
                reTrmFormulas = reTrmFormulas,
                reTrmAddArgs  = reTrmAddArgs,
                reTrmClasses  = reTrmClasses))
}

reParen <- function(reTrm) paste("(", deparse(reTrm), ")", sep = "", collapse = "")

##' @rdname splitForm
##' @param splitFormula results of \code{splitForm}
##' @export
reForm <- function(splitFormula) {
    characterPieces <- c(list(deparse(splitFormula$fixedFormula)),
                         lapply(splitFormula$reTrmFormulas, reParen))
    as.formula(do.call(paste, c(characterPieces, list(sep = " + "))))
}


##' @rdname splitForm
##' @export
noSpecials <- function(term) {
    nospec <- noSpecials_(term)
    if (is(term,"formula") && length(term)==3 && is.symbol(nospec)) {
        ## called with two-sided RE-only formula:
        ##    construct response~1 formula
        nospec <- reformulate("1", response = deparse(nospec))
    }
    return(nospec)
}

noSpecials_ <- function(term) {
    if (!anySpecial(term)) return(term)
    if (isSpecial(term)) return(NULL)
    nb2 <- noSpecials(term[[2]])
    nb3 <- noSpecials(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

isSpecial <- function(term) {
    if(is.call(term)) {
        for(cls in findReTrmClasses()) {
            if(term[[1]] == cls) return(TRUE)
        }
    }
    FALSE
}

isAnyArgSpecial <- function(term) {
    for(i in seq_along(term)) {
        if(isSpecial(term[[i]])) return(TRUE)
    }
    FALSE
}

anySpecial <- function(term) {
    any(findReTrmClasses() %in% all.names(term))
}

##' Make strucGlmer object
##'
##' @param opt,parsedForm,dfun,mc See \code{\link{strucGlmer-class}}
##' @return An object of \code{\link{strucGlmer-class}}
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

