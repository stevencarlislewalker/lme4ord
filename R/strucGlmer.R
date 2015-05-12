##' Structured glmer
##'
##' Parses a mixed model formula with special structures (with
##' \code{\link{strucParseFormula}}), constructs a generalized linear
##' mixed model deviance function (with
##' \code{\link{mkGeneralGlmerDevfun}}), optimizes this deviance
##' function (with \code{\link{bobyqa}}), and returns the results as a
##' \code{strucGlmer} object.
##'
##' @param formula extended mixed model formula
##' @param data data frame
##' @param addArgs list of additional arguments to pass to
##' @param optVerb verbose
##' \code{\link{setReTrm}} methods
##' @param ... further arguments to \code{\link{mkGeneralGlmerDevfun}}
##' @export
strucGlmer <- function(formula, data, addArgs = list(), optVerb = 0L,
                         weights, offset, etastart,
                         ...) {

    cat("\nConstructing vectors and matrices...\n")
    pForm <- strucParseFormula(formula, data, addArgs, ...)

    cat("\nConstructing deviance function...\n")
    dfun <- mkGeneralGlmerDevfun(y = pForm$response,
                                 X = pForm$fixed,
                                      Zt = as(pForm$Zt,      "dgCMatrix"),
                                 Lambdat = as(pForm$Lambdat, "dgCMatrix"),
                                 weights = weights,
                                 offset = offset,
                                 etastart = etastart,
                                 initPars = pForm$initPars,
                                 parInds = pForm$parInds,
                                 mapToCovFact = mkSparseTrans(pForm$Lambdat),
                                 mapToModMat = mkSparseTrans(pForm$Zt),
                                 devfunEnv = pForm$devfunEnv,
                                 ...)

    cat("\nInitializing deviance function...\n")
    dfun(pForm$initPars)

    cat("\nOptimizing deviance function...\n")
    opt <- minqa:::bobyqa(pForm$initPars, dfun, lower = pForm$lower,
                          control =
                          list(iprint = optVerb,
                               rhobeg = 0.0002,
                               rhoend = 2e-7))

    cat("\nPreparing output...\n")
    names(opt$par) <- names(pForm$initPars)

    ans <- list(opt = opt, parsedForm = pForm, dfun = dfun)
    class(ans) <- "strucGlmer"
    
    try(ansRand <- mapply(setInit, pForm$random,
                          covarPerTerm(ans),
                          loadsPerTerm(ans),
                          SIMPLIFY = FALSE), silent = TRUE)
    if(inherits(ansRand, "try-error")) {
        warning("couldn't reset initial values to estimated values,\n",
                "so model printout is probably misleading.\n",
                "perhaps check the object itself.")
    } else {
        ans$pForm$random <- ansRand
    }
    return(ans)
}

##' @rdname strucGlmer
##' @export
print.strucGlmer <- function(x, ...) {
    cat ("Structured generalized linear mixed model\n")
    cat ("=========================================\n")
    cat ("\nFixed effects\n")
    cat (  "-------------\n")
    print(fixef(x))
    cat ("\nRandom effects\n")
    cat (  "--------------\n")
    lapply(x$parsedForm$random, printReTrm)
}

##' @param type character string giving the type of parameter
##' (e.g. \code{"fixef", "covar"})
##' @rdname strucGlmer
##' @export
getStrucGlmerPar <- function(object, type, ...) {
    parInds <- environment(object$dfun)$parInds
    optPar <- object$opt$par
    optPar[unlist(parInds[type])]
}

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

##' @rdname strucGlmer
##' @export
loads.strucGlmer <- function(object, ...) {
    getStrucGlmerPar(object, "loads")
}

subRagByLens <- function(x, lens) {
    if(is.null(names(lens))) names(lens) <- seq_along(lens)
    split(x, rep(names(lens), lens))
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
##' @rdname strucParseFormula
##' @export
strucParseFormula <- function(formula, data, addArgs = list(), reTrmsList = NULL, ...) {
                                        # get and construct basic
                                        # information: (1) list of
                                        # formulas, (2) data, and (3)
                                        # the environment that will
                                        # eventually become the
                                        # environment of the deviance
                                        # function
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

    ## FIXME: more flexibility for lower, and maybe add an upper?
    lower <- c(ifelse(init$covar, 0, -Inf),
               rep(-Inf, length(initPars) - length(init$covar)))
    names(lower) <- names(initPars)

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
                lower = lower,
                devfunEnv = devfunEnv))
}


##' @param parsedForm result of \code{strucParseFormula}
##' @param family family object
##' @param weights optional weights
##' @rdname strucParseFormula
##' @export
simStrucParsedForm <- function(parsedForm, family = binomial,
                                 weights, ...) {
    if(missing(weights)) weights <- rep(1, length(parsedForm$response))
    with(parsedForm, {
        reMM <- t(as(Zt, "dgCMatrix")) %*% t(as(Lambdat, "dgCMatrix"))
        feMM <- fixed
        fe <- as.numeric(feMM %*% initPars[parInds$fixef])
        re <- as.numeric(reMM %*% rnorm(ncol(reMM)))
        simFun <- simfunList[[family()$family]]
        return(simFun(weights, length(weights), family()$linkinv(fe + re)))
    })
}


##' Make general random effects terms
##'
##' @param reStructList list of random effects structure functions
##' @param parsedForm list parsed random effects terms
##' @param ... potential extra parameters
##' @rdname mkReTrms
##' @export
mkGeneralReTrms <- function(reStructList, parsedForm, ...) {
    trmList <- listTranspose(mapply(do.call,
                                    reStructList,
                                    parsedForm$random,
                                    SIMPLIFY = FALSE))
         ZtBind <- as.repSparse(sort(.bind(trmList$Zt,      "row" )))
    LambdatBind <- as.repSparse(sort(.bind(trmList$Lambdat, "diag")))
    return(list(Zt = as(ZtBind, "dgCMatrix"),
                Lambdat = as(LambdatBind, "dgCMatrix"),
                ZtTrans = mkSparseTrans(ZtBind),
                LambdatTrans = mkSparseTrans(LambdatBind)))
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

                                        # find additional arguments
    reTrmAddArgs <- lapply(formSplits, "[", -2)[!(formSplitID == "(")]
                                        # change call name
    reTrmAddArgs <- lapply(reTrmAddArgs, "[[<-", 1, as.name("c"))
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
    #if (isAnyArgSpecial(term)) return(NULL)
    #if (length(term) == 2) {
    #    term[[2]] <- noSpecials(term[[2]])
    #    return(term)
    #}
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


##' Version of the recursive nobars function from lme4
##'
##' @param term term
##' @export
nobarsWithSpecials <- function (term) {
    specials <- findReTrmClasses()
    if(length(term) > 2L) {
        if(any(as.character(term[[1]]) == specials)) {
            term <- term[1:2]
        }
    }
    ant <- all.names(term)
    if ((!any(c("|", "||") %in% ant)) && (!any(ant %in% specials)))
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (is.call(term) && term[[1]] == as.name("||")) 
        return(NULL)
    #if (as.character(term[[1]]) %in% specials) 
    #    term <- term[[2]]
    if (length(term) == 1)
        return(term)
    if (length(term) == 2) {
        nb <- nobarsWithSpecials(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobarsWithSpecials(term[[2]])
    nb3 <- nobarsWithSpecials(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}


## ##' Get components from a structured glmer model
## ##'
## ##' @param object a \code{\link{strucGlmer}} object
## ##' @param name a character vector naming the component
## ##' @export
## getSG <- function(object,
##                   name =
##                   c("y", "X", "Zt", "Lambdat",
##                     "")) {
## }


