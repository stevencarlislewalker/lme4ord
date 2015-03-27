generalGlmer <- function(formula, data, ...) {
    generalParseFormula(formula, data, ...)
}

##' Parse a mixed model formula
##'
##' @param formula mixed model formula
##' @param data an object coercible to data frame
##' @param ... additional parameters to \code{\link{as.data.frame}}
##' @export
generalParseFormula <- function(formula, data, ...) {
    formula <- splitForm(formula)
    data <- as.data.frame(data, ...)

    with(formula, mapply(mkReTrmStruct,
                         reTrmFormulas,
                         reTrmClasses,
                         SIMPLIFY = FALSE))
    
                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

    return(list(response = model.response(model.frame(nobars(formula), data)),
                fixed    = model.matrix(nobars(formula), data),
                random   = reTrmsList))
}

##' Make general random effects terms
##'
##' @param reStructList list of random effects structure functions
##' @param parsedForm list parsed random effects terms
##' @param ... potential extra parameters
##' @rdname mkReTrms
##' @export
mkGeneralReTrms <- function(reStructList, parsedForm, ...) {
    trmList <- mapply(do.call, reStructList, parsedForm$random, SIMPLIFY = FALSE)
    trmList <- listTranspose(trmList)
         ZtBind <- as.repSparse(sort(.bind(trmList$Zt,      "row" )))
    LambdatBind <- as.repSparse(sort(.bind(trmList$Lambdat, "diag")))
    return(list(Zt = as(ZtBind, "dgCMatrix"),
                Lambdat = as(LambdatBind, "dgCMatrix"),
                ZtTrans = mkSparseTrans(ZtBind),
                LambdatTrans = mkSparseTrans(LambdatBind)))
}

##' Make transformation function for column-compressed sparse matrix
##'
##' @param object repeated sparse matrix object
##' @rdname mkReTrms
##' @export
mkSparseTrans <- function(object) {
    local({
        trans <- object$trans
        inds <- object$valInds
        function(matPars) trans(matPars)[inds]
    })
}

##' Find classes with a \code{getReTrm} method
##'
##' @param formula generalized mixed model formula.  if \code{NULL}
##' (the default) \code{findReTrmClasses} returns classes available
##' (on the search path).
##' @export
findReTrmClasses <- function(formula = NULL) {
    if(is.null(formula)) {
        return(as.character(sub("getReTrm.", "", methods("getReTrm"))))
    }
    classInds <- attr(terms(formula, specials = findReTrmClasses()), "specials")
    unlist(mapply(rep, names(classInds),
                  lapply(classInds, length),
                  SIMPLIFY = FALSE))[unlist(classInds)]
}


##' Split a formula
##'
##' @param formula Generalized mixed model formula
##' @export
splitForm <- function(formula) {

    specials <- findReTrmClasses()

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
                                        # without a getReTrm method
    formSplitID <- sapply(lapply(formSplits, "[[", 1), as.character)
                                        # warn about terms without a
                                        # getReTrm method
    badTrms <- formSplitID == "|"
    if(any(badTrms)) {
        warning(paste("can't find getReTrm method(s) for term number(s)",
                      paste(which(badTrms), collapse = ", "),
                      "\ntreating those terms as unstructured"))
        formSplitID[badTrms] <- "("
        fixBadTrm <- function(formSplit) {
            as.formula(paste(c("~(", as.character(formSplit)[c(2, 1, 3)], ")"),
                             collapse = " "))[[2]]
        }
        formSplits[badTrms] <- lapply(formSplits[badTrms], fixBadTrm)
    }
                                        # standard RE terms
    formSplitStan <- formSplits[formSplitID == "("]
                                        # structured RE terms
    formSplitSpec <- formSplits[!(formSplitID == "(")]

    if(length(formSplitSpec) == 0) stop(
                 "no special covariance structures. ",
                 "please use lmer or glmer")


    fixedFormula <- formula(paste(formula[[2]], "~", as.character(nobars(formula))[[3]]))
    reTrmFormulas <- c(lapply(formSplitStan, "[[", 2),
                       lapply(formSplitSpec, "[[", 2))
    reTrmClasses <- c(rep("reTrmFull", length(formSplitStan)),
                      sapply(lapply(formSplitSpec, "[[", 1), as.character))
    
    return(list(fixedFormula = fixedFormula,
                reTrmFormulas = reTrmFormulas,
                reTrmClasses = reTrmClasses))
}

##' Make random effects term structure
##'
##' @param formula random effects term formula (without special function)
##' @param class random effects term class
##' @param data a data frame
##' @export
mkReTrmStruct <- function(formula, class, data) {
    ans <- getModMatAndGrpFac(formula, data)
    names(ans) <- ans[["grpName"]] # FIXME: include struct class name?
    ans$formula <- formula
    structure(ans, class = class)
}
