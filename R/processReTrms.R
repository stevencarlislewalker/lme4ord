##' Set model matrix slice and relative covariance factor block for a
##' random effects term
##' 
##' @param object a \code{reTrmStruct} object
##' @param addArgs named list of additional arguments
##' @param devfunEnv optional environment of the deviance function
##' @rdname setReTrm
##' @seealso \code{\link{mkReTrmStructs}} for construction of these objects
##' @return \code{object} with the following additional elements:
##' 
##' \item{Zt}{A \code{\link{repSparse}} object describing the slice of
##' the model matrix associated with the random effects term.}
##'
##' \item{Lambdat}{A \code{\link{repSparse}} object describing the
##' block of the relative covariance factor associated with the random
##' effects term.}
##' @export
setReTrm <- function(object, addArgs, devfunEnv = NULL) {
    UseMethod("setReTrm")
}

##' @rdname setReTrm
##' @export
setReTrm.default <- function(object, addArgs, devfunEnv = NULL) {
    Zt <- kr(t(as.repSparse(object$modMat)), as.repSparse(object$grpFac))
    nc <- ncol(object$modMat)
    nl <- nlevels(object$grpFac)
    template <- repSparseTri(rep(1, nc), rep(0, choose(nc, 2)), FALSE)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = rep(template, nl, type = "diag"))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.identity <- function(object, addArgs, devfunEnv = NULL) {
    Zt <- kr(t(as.repSparse(object$modMat)), as.repSparse(object$grpFac))
    nl <- nlevels(object$grpFac)
    nc <- ncol(object$modMat)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nl * nc))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.flexvar <- function(object, addArgs, devfunEnv = NULL) {
    n <- nrow(object$modMat)
    addArgs <- eval(object$addArgs, addArgs)
    Zt <- repSparseIdent(n)
    ##mkFlexObsLevelTrans(
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(n))),
                     class = class(object)))
##    addArgs <- eval(object$addArgs, addArgs)
##    Zt <- kr(t(as.repSparse(object$modMat)), as.repSparse(object$grpFac))  
}

##' @rdname setReTrm
##' @export
setReTrm.edge <- function(object, addArgs, devfunEnv = NULL) {
    addArgs <- eval(object$addArgs, addArgs)
    Jedge <- as(edgeTipIndicator(addArgs$phy), "sparseMatrix")
    Jspp <- as(object$grpFac, "sparseMatrix")
    Jt <- as.repSparse(Jedge %*% Jspp)
    Zt <- kr(t(as.repSparse(object$modMat)), Jt)
    nl <- nrow(Jedge)
    nc <- ncol(object$modMat)
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = repSparseIdent(nl * nc))),
                     class = class(object)))
}

##' @rdname setReTrm
##' @export
setReTrm.cooccur <- function(object, addArgs, devfunEnv = NULL) {
    addArgs <- eval(object$addArgs, addArgs)
    Jt <- as.repSparse(as(object$grpFac, "sparseMatrix"))
    Zt <- kr(t(as.repSparse(object$modMat)), Jt)
    nCovPars <- choose(ncol(object$modMat), 2)
    Tt <- t(repSparseCorMatChol(rep(0, nCovPars)))
    return(structure(c(object,
                       list(Zt = resetTransConst(Zt),
                            Lambdat = rep(Tt, nrow(Jt), type = "diag")),
                       class = class(object))))
}

##' Update a random effects term structure with new parameters
##'
##' @param object a \code{reTrmStruct} object
##' @param newCovar new covariance parameters
##' @param newLoads new loadings parameters
##' @export
update.reTrmStruct <- function(object, newCovar, newLoads) {
    setInit(object$Lambdat, newCovar)
    setInit(object$Zt,      newLoads)
    return(object)
}

##' Print random effects term
##'
##' @param object \code{\link{repSparse}} object
##' @param forSummary print for \code{\link{summary}} instead of
##' \code{\link{print}}?
##' @export
printReTrm <- function(object, forSummary = FALSE, ...) {
    UseMethod("printReTrm")
}

.printPars <- function(description = "parameters: ", value) {
    if((length(value) > 0L) && (!is.na(value)) && (!is.null(value))) {
        cat(description, value)
    }
}

##' @rdname printReTrm
##' @export
printReTrm.default <- function(object, forSummary = FALSE, ...) {
    cat(paste("A", class(object)[1], "random effects structure\n"))
    .printPars("  covariance parameters: ", getInit(object$Lambdat))
    .printPars("  loadings parameters:   ", getInit(object$Zt))
    cat("\n")
}


##' Simulate additional arguments
##'
##' @param object a \code{reTrmStruct} object
##' @param ... dots
##' @rdname simAddArgs
##' @export
simAddArgs <- function(object, ...) {
    UseMethod("simAddArgs")
}

##' @rdname simAddArgs
##' @export
simAddArgsList <- function(object, ...) {
    l... <- list(...)
    unlist(do.call(lapply, c(list(unname(object), simAddArgs), l...)), FALSE)
}

##' @rdname simAddArgs
##' @export
simAddArgs.default <- function(object, ...) list()

##' @param rtreeArgs arguments for \code{\link{rtree}}
##' @param compute.brlenArgs arguments for \code{\link{compute.brlen}}
##' @rdname simAddArgs
##' @export
simAddArgs.edge <- function(object, rtreeArgs = list(),
                            compute.brlenArgs = list(), ...) {
    namePhy <- as.character(object$addArgs$phy)
    phy <- do.call(rtree, c(list(nlevels(object$grpFac)), rtreeArgs))
    phy <- do.call(compute.brlen, c(list(phy), compute.brlenArgs))
    phy$tip.label <- unique(as.character(object$grpFac))
    return(setNames(list(phy), namePhy))
}


##' @param formula generalized mixed model formula.  if \code{NULL}
##' (the default) \code{findReTrmClasses} returns classes available
##' (on the search path).
##' @rdname setReTrm
##' @export
findReTrmClasses <- function(formula = NULL) {
    if(is.null(formula)) {
        return(as.character(sub("setReTrm.", "", methods("setReTrm"))))
    }
    ## intersect(all.names(formula), findReTrmClasses())
    classInds <- attr(terms(formula, specials = findReTrmClasses()), "specials")
    names(unlist(classInds))
    ## unlist(mapply(rep, names(classInds),
    ##               lapply(classInds, length),
    ##               SIMPLIFY = FALSE))[unlist(classInds)]
}
