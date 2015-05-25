##' Set model matrix slice and relative covariance factor block for a
##' random effects term
##' 
##' @param object a \code{reTrmStruct} object
##' @param addArgsList a list of named quantities within which
##' \code{addArgsExpr} is evaluated
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
##'
##' \item{lowerLoads,upperLoads,lowerCovar,upperCovar}{Lower/upper
##' bounds on loading parameters (possible parameters for \code{Zt})
##' and/or covariance parameters (possible parameters for \code{Lambdat})}
##'
##' @export
setReTrm <- function(object, addArgsList,
                     devfunEnv = NULL) {
    UseMethod("setReTrm")
}

##' @rdname setReTrm
##' @export
setReTrm.default <- function(object, addArgsList,
                             devfunEnv = NULL) {
    
                                        # transposed model matrix (or
                                        # loadings matrix) -- Zt rows
                                        # associated with this term
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)),
                               as.repSparse(object$grpFac)))

                                        # covariance factor -- block
                                        # of Lambdat associated with
                                        # this term
    nc <-    ncol(object$modMat)
    nl <- nlevels(object$grpFac)
    templateBlock <- repSparseTri(   diagVals = rep(1,        nc    ),
                                  offDiagVals = rep(0, choose(nc, 2)),
                                  low = FALSE)
    Lambdat <- rep(templateBlock, nl, type = "diag")
    attr(Lambdat, "templateBlock") <- templateBlock ## for VarCorr

                                        # package up object
                                        # (implicitly sets lower and
                                        # upper bounds, but these can
                                        # be explicitly set with
                                        # packReTrm)
    packReTrm(object, Zt, Lambdat)
}

##' @rdname setReTrm
##' @export
setReTrm.factAnal <- function(object, addArgsList,
                              devfunEnv = NULL) {
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    nl <- nlevels(grpFac <- object$grpFac)
    nc <- addArgs$nAxes
    set.seed(addArgs$seed)
    modMat <- subset(rRepSparse(nl,  nc,
                                nl * nc,
                                nl * nc), as.numeric(grpFac))
    
    Zt <- kr(t(modMat), resetTransConst(as.repSparse(addArgs$obsFac)))
    Lambdat <- resetTransConst(repSparseIdent(nrow(Zt)))
    
    lowerLoads <- rep(-Inf, length(getInit(Zt)))
    upperLoads <- rep( Inf, length(getInit(Zt)))

    packReTrm(object, Zt, Lambdat,
              lowerLoads = lowerLoads,
              upperLoads = upperLoads)
}


##' @rdname setReTrm
##' @export
setReTrm.identity <- function(object, addArgsList, devfunEnv = NULL) {

                                        # Zt
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)),
                             as.repSparse(object$grpFac)))

                                        # Lambdat
    nl <- nlevels(object$grpFac)
    nc <-    ncol(object$modMat)
    Lambdat <- repSparseIdent(nl * nc)

                                        # pack
    packReTrm(object, Zt, Lambdat)
}

##' @rdname setReTrm
##' @export
setReTrm.flexvar <- function(object, addArgsList, devfunEnv = NULL) {
                                        # get additional arguments and
                                        # sample size
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    n <- nrow(object$modMat)

                                        # Zt
    Zt <- resetTransConst(repSparseIdent(n))

                                        # Lambdat
    inds <- seq_len(n); baselineVars <- rep(1, n)
    if(is.null(init <- addArgs$init)) init <- rep(0, addArgs$nBasis)
    Lambdat       <- repSparseDiag  (baselineVars, inds)
    Lambdat$trans <- mkFlexDiagTrans(init, baselineVars, devfunEnv)

                                        # pack
    packReTrm(object, Zt, Lambdat,
              lowerCovar = rep(-Inf, length(init)))
}


##' @rdname setReTrm
##' @export
setReTrm.expcutdist <- function(object, addArgsList, devfunEnv = NULL) {
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    distCutoff <- addArgs$distCutoff
    minCov <- addArgs$minCov
    distMat <- edgeMat <- addArgs$distMat
    
    edgeMat[] <- 1 * (distMat[] < distCutoff)
    sparseMat <- as(as.matrix(distMat), "TsparseMatrix")
    inds <- (sparseMat@x < distCutoff) & (sparseMat@i > sparseMat@j)
    Jedge <- sparseMatrix(i = rep(seq_along(which(inds)), 2),
                          j = c(sparseMat@i[inds] + 1, sparseMat@j[inds] + 1),
                          x = rep(1, 2 * sum(inds)))
    colnames(Jedge) <- attr(distMat, "Labels")
    edgeDists <- sparseMat@x[inds]
    Jgrp <- as(object$grpFac, "sparseMatrix")
    Jt <- as.repSparse(Jedge %*% Jgrp)
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)), Jt))
    transFun <- local({
        init <- c(0.01)
        edgeDists <- edgeDists
        distCutoff <- distCutoff
        minCov <- minCov
        function(matPars) {
            q1 <- (minCov - 1)/(exp(-(matPars[1]) * distCutoff) - 1)
            q2 <- 1 - q1
            (q2 + q1 * exp(-(matPars[1]) * edgeDists))
        }
    })
    Lambdat <- repSparseDiag(transFun(1))
    Lambdat$trans <- transFun
    packReTrm(object, Zt, update(Lambdat))
}

##' @rdname setReTrm
##' @export
setReTrm.edge <- function(object, addArgsList, devfunEnv = NULL) {
                                        # get additional arguments
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)

                                        # Zt
    Jedge <- as(edgeTipIndicator(addArgs$phy), "sparseMatrix")
    Jspp <- as(object$grpFac, "sparseMatrix")
    Jt <- resetTransConst(as.repSparse(Jedge %*% Jspp))
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)), Jt))

                                        # Lambdat
    nl <- nrow(Jedge)
    nc <- ncol(object$modMat)
    Lambdat <- repSparseIdent(nl * nc)
    
                                        # pack
    packReTrm(object, Zt, Lambdat)
}

##' @rdname setReTrm
##' @export
setReTrm.cooccur <- function(object, addArgsList, devfunEnv = NULL) {

                                        # get additional arguments
    addArgs <- eval(object$addArgs, addArgsList)

                                        # Zt
    Jt <- as.repSparse(as(object$grpFac, "sparseMatrix"))
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)), Jt))

                                        # Lambdat
    nCovPars <- choose(ncol(object$modMat), 2)
    Tt <- t(repSparseCorMatChol(rep(0, nCovPars)))
    Lambdat <- rep(Tt, nrow(Jt), type = "diag")
    
                                        # pack
    packReTrm(object, Zt, Lambdat)
}

##' @rdname setReTrm
##' @export
setReTrm.varIdent <- function(object, addArgsList, devfunEnv = NULL) {
                                        # transposed model matrix (or
                                        # loadings matrix) -- Zt rows
                                        # associated with this term
    n <- length(object$grpFac)
    Zt <- resetTransConst(repSparseIdent(n))
    
                                        # covariance factor -- block
                                        # of Lambdat associated with
                                        # this term

    nl <- nlevels(object$grpFac)
    init <- rep(1, nl)
    Lambdat <- repSparseVarWithCovariate(init, 
                                         grpFac = object$grpFac,,
                                         mkTransFunc = mkVarIdentTrans)
    
                                        # package up object
    packReTrm(object, Zt, Lambdat,
              lowerCovar = rep(0, nl))
}

##' @rdname setReTrm
##' @export
setReTrm.multiVarExp <- function(object, addArgsList, devfunEnv = NULL) {
                                        # transposed model matrix (or
                                        # loadings matrix) -- Zt rows
                                        # associated with this term
    n <- nrow(object$modMat)
    Zt <- resetTransConst(repSparseIdent(n))
    
                                        # covariance factor -- block
                                        # of Lambdat associated with
                                        # this term

    if(is.na(object$grpFac)) {
        nl <- 1
    } else {
        nl <- nlevels(object$grpFac)
    }
    nc <- ncol(object$modMat)
    init <- rep(0, nl * nc)
    Lambdat <- repSparseDiag(rep(1, n), 1:n)
    Lambdat$trans <- mkMultiVarExpTrans(init, object$modMat, object$grpFac)
    
                                        # package up object
    packReTrm(object, Zt, Lambdat,
              lowerCovar = rep(-Inf, length(init)))
}

##' @rdname setReTrm
##' @export
setReTrm.obslev <- function(object, addArgsList,
                             devfunEnv = NULL) {

    n <- nrow(object$modMat)
    Zt <- resetTransConst(repSparseIdent(n))
    Lambdat <- repSparseIdent(n)
    packReTrm(object, Zt, Lambdat)
}


##' @param Zt transposed model matrix
##' @param Lambdat relative covariance factor
##' @param lowerLoads,upperLoads,lowerCovar,upperCovar lower and upper
##' bounds on possible loadings parameters (for \code{Zt}) and
##' possible covariance parameters (for \code{Lambdat})
##' @rdname setReTrm
##' @export
packReTrm <- function(object, Zt, Lambdat,
                      lowerLoads, upperLoads,
                      lowerCovar, upperCovar) {
    if(missing(lowerLoads)) lowerLoads <- setLowerDefault(getInit(Zt))
    if(missing(upperLoads)) upperLoads <- setUpperDefault(getInit(Zt))
    if(missing(lowerCovar)) lowerCovar <- setLowerDefault(getInit(Lambdat))
    if(missing(upperCovar)) upperCovar <- setUpperDefault(getInit(Lambdat))
    structure(c(object,
                list(Zt = Zt, Lambdat = Lambdat,
                     lowerLoads = lowerLoads,
                     upperLoads = upperLoads,
                     lowerCovar = lowerCovar,
                     upperCovar = upperCovar)),
              class = class(object))
}

##' @export
VarCorr.reTrmStruct <- function(x, sigma = 1, rdig = 3) {
    ## default method
    if(!is.null(templateBlock <- attr(x$Lambdat, "templateBlock"))) {
        ans <- update(templateBlock, getInit(x$Lambdat))
        vc <- as.matrix(Matrix::crossprod(as.matrix(ans, sparse = TRUE)))
        rownames(vc) <- colnames(vc) <- colnames(x$modMat)
    } else {
        vc <- as.matrix(Matrix::crossprod(as.matrix(x$Lambdat, sparse = TRUE)))
    }
    attr(vc, "stddev") <- sqrt(diag(as.matrix(vc)))
    attr(vc, "correlation") <- cov2cor(vc)
    attr(vc, "useSc") <- FALSE
    return(vc)
}

##' Update a random effects term structure with new parameters
##'
##' @param object a \code{reTrmStruct} object
##' @param newCovar new covariance parameters
##' @param newLoads new loadings parameters
##' @param ... potential additional arguments (ignored by the default
##' method)
##' @rdname update.reTrmStruct
##' @export
update.reTrmStruct <- function(object, newCovar, newLoads, ...) {
    object$Lambdat <- update(object$Lambdat, newCovar)
    object$Zt      <- update(object$Zt,      newLoads)
    return(object)
}

##' @rdname update.reTrmStruct
##' @export
update.flexvar <- function(object, newCovar, newLoads, ...) {

    ## This special method is to put information in the environment of
    ## the transformation function that is currently only in the
    ## environment of the deviance function.  _In general_ such a
    ## method is required whenever the transformation function depends
    ## on other random effects terms.  Please see `?assignWith` for a
    ## useful technique in this regard
    
    object <- update.reTrmStruct(object, newCovar, newLoads)
    transEnv <- environment(object$Lambdat$trans)
    assignWith(expr  = indsForClass("flexvar", reTrmClasses, nRePerTrm),
               name  = "indsObsLevel",
               data  = transEnv$devfunEnv,
               envir = transEnv)
    return(object)
}

##' Get defaults choices for lower and/or upper bound of a model
##' parameter
##'
##' @param init initial parameter value
##' @param ... additional arguments not currently used
##' @rdname setLowerUpperDefault
##' @export
setLowerDefault <- function(init, ...) {
    if(missing(init)) return(NULL)
    if(length(init) == 0L) return(NULL)
    ifelse(as.logical(init), 0, -Inf)
}


##' @rdname setLowerUpperDefault
##' @export
setUpperDefault <- function(init, ...) {
    if(missing(init)) return(NULL)
    if(length(init) == 0L) return(NULL)
    rep(Inf, length(init))
}


indsForClass <- function(reTrmClass, reTrmClasses, nValuesPerTrm) {
    ## FIXME: not efficient -- computes indices for all classes first
    whichClass <- which(reTrmClasses == reTrmClass)
    starts <- c(1L, 1L + cumsum(nValuesPerTrm)[-length(nValuesPerTrm)])
    ends <- starts + nValuesPerTrm - 1L
    unlist(mapply(":", starts[whichClass], ends[whichClass], SIMPLIFY = FALSE))
}

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

.printPars <- function(description = "parameters: ", value) {
    if((length(value) > 0L) && (!is.na(value)) && (!is.null(value))) {
        cat(description, value, "\n")
    }
}

.printVC <- function(description = "variance-correlation: ", value) {
    if((nrow(value[[1]]) < 6) && (!is.null(rownames(value[[1]])))) {
        cat(description, "\n")
        print(lme4:::formatVC(value), quote = FALSE, digits = 3)
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

## FIXME: write specific printReTrm methods for different classes


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
    classInds <- attr(terms(formula, specials = findReTrmClasses()), "specials")
    names(unlist(classInds))
}

##' @param addArgsExpr a list of expressions for evaluating within
##' \code{addArgsList}
##' @rdname setReTrm
##' @export
getAddArgs <- function(addArgsExpr, addArgsList) {
    setNames(lapply(addArgsExpr, eval, envir = addArgsList),
              names(addArgsExpr))
}

##' Simulate from a random effects term (experimental)
##'
##' @param object \code{reTrmStruct} object
##' @export
simReTrm <- function(object) {
    with(object, {
        as.numeric(rnorm(nrow(Lambdat)) %*% as.matrix(Lambdat) %*% as.matrix(Zt))
    })
}
