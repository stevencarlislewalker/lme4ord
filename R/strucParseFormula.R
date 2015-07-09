## ----------------------------------------------------------------------
## Formula parsing
##
## used by:  strucGlmer.R
## uses:  
## ----------------------------------------------------------------------


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
##' \item{mapToModMat}{function taking the \code{loads} parameters and
##' returning the values of the non-zero elements of \code{Zt}
##' (i.e. the \code{x} slot of \code{Zt})}
##' \item{mapToCovFact}{function taking the \code{covar} parameters
##' and returning the values of the non-zero elements of
##' \code{Lambdat} (i.e. the \code{x} slot of \code{Lambdat})}
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
    ZtBind      <- .bind(     ZtList,  "row")
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
                mapToModMat = mkSparseTrans(Zt),
                mapToCovFact = mkSparseTrans(Lambdat),
                initPars = initPars, parInds = parInds,
                lower = lower, upper = upper,
                devfunEnv = devfunEnv,
                formula = formula))
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
                 "please use lmer or glmer, ",
                 "or use findReTrmClasses() for available structures.")


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


##' @param term language object
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

