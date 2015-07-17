## ----------------------------------------------------------------------
## Random-effects structures
##
## used by:  strucParseFormula.R, outputModule.R
## uses:  repSparse.R
## ----------------------------------------------------------------------


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
    if(inherits(splitFormula, "formula")) {
        splitFormula <- splitForm(splitFormula)
    }
    reTrmsList <- lapply(splitFormula$reTrmFormulas,
                         getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- paste(sapply(reTrmsList, "[[", "grpName"),
                               splitFormula$reTrmClasses, sep = ".")
    nUnStr <- sum(splitFormula$reTrmClasses == "unstruc")
    for(i in seq_along(reTrmsList)) {
        clsi <- splitFormula$reTrmClasses[[i]]
        if(clsi != "unstruc") {
            reTrmsList[[i]]$addArgs <- splitFormula$reTrmAddArgs[[i - nUnStr]]
        }
        class(reTrmsList[[i]]) <- c(clsi, "reTrmStruct")
    }
    return(reTrmsList)
}

##' Get model matrix and grouping factor
##'
##' This is kind of like \code{\link{model.matrix}} but for random
##' effects terms.
##' 
##' @param bar random effect language object (e.g. \code{x | g})
##' @param fr model frame
##' @return list with model matrix and grouping factor
##' @export
getModMatAndGrpFac <- function(bar, fr) {
    ## based on mkBlist

    noGrpFac <- is.null(lme4::findbars(bar))

    if(!noGrpFac) {
        linFormLang <- bar[[2]] # language object specifying linear model
            grpLang <- bar[[3]] # language object specifying grouping factor

        fr <- lme4::factorize(bar, fr)
        nm <- deparse(grpLang)
        ## try to evaluate grouping factor within model frame ...
        if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac),
                                                   list(fac = grpLang)), fr),
                                   error = function(e) NULL)))
            stop("couldn't evaluate grouping factor ",
                 nm, " within model frame:",
                 " try adding grouping factor to data ",
                 "frame explicitly if possible", call. = FALSE)
        if (all(is.na(ff)))
            stop("Invalid grouping factor specification, ",
                 nm, call. = FALSE)
    } else { # noGrpFac
        linFormLang <- bar
        ff <- nm <- NA
    }
    mm <- model.matrix(eval(substitute( ~ foo, list(foo = linFormLang))), fr)
    return(list(modMat = mm, grpFac = ff, grpName = nm))
}


##' Set model matrix slice and relative covariance factor block for a
##' random effects term
##' 
##' @param object a \code{reTrmStruct} object
##' @param addArgsList a list of named quantities within which
##' \code{addArgsExpr} is evaluated
##' @param auxEnv an optional auxilliary environment containing
##' objects that are possibly required for the setting of a random
##' effects term structure (currently this is almost always the
##' environment of a call to \code{\link{strucParseFormula}}, unless
##' \code{setReTrm} is called directly)
##' @param devfunEnv optional environment of the deviance function
##' @rdname setReTrm
##' @note Almost never called directly, but instead called within
##' \code{\link{strucParseFormula}}.
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
##' @section Mixed model formula usage:
##' The default method is identical to a standard \code{lme4} random effects term.
##' \code{(linForm | grpFac)}
##' \describe{
##' \item{linForm}{linear model formula}
##' \item{grpFac}{grouping factor}
##' }
##'
##' @seealso getReTrm
##' @family setReTrm
##' @export
setReTrm <- function(object, addArgsList,
                     auxEnv = NULL, devfunEnv = NULL) {
    UseMethod("setReTrm")
}

##' @rdname setReTrm
##' @export
setReTrm.default <- function(object, addArgsList,
                             auxEnv = NULL, devfunEnv = NULL) {
    
                                        # transposed model matrix (or
                                        # loadings matrix) -- Zt rows
                                        # associated with this term
    Zt <- resetTransConst(kr(as.repSparse(object$grpFac),
                           t(as.repSparse(object$modMat))))

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
    packReTrm(object, Zt, Lambdat, devfunEnv = devfunEnv)
}

##' lme4 emulator
##' 
##' @export
##' @template setReTrm
##' @templateVar cls lme4
##' @templateVar form \code{lme4(linForm | grpFac)}
##' @templateVar arg "grpFac"
##' @templateVar desc "grouping factor"
setReTrm.lme4 <- function(object, addArgsList,
                          auxEnv = NULL, devfunEnv = NULL) {
    setReTrm.default(object, addArgsist, devfunEnv)
}

##' Random effects term for factor analysis
##'
##' @export
##' @template setReTrm
##' @templateVar cls factAnal
##' @templateVar form \code{factAnal(0 + obsFac | varFac, nAxes, seed)}
##' @templateVar arg c("obsFac", "varFac", "nAxes", "seed")
##' @templateVar desc c("grouping factor for multivariate observations", "grouping factor for variables", "number of axes", "random seed for initial axis values")
setReTrm.factAnal <- function(object, addArgsList,
                              auxEnv = NULL, devfunEnv = NULL) {
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)

    nVar <- nlevels(varFac <- object$grpFac) # number of variables
    nObs <- ncol(obsMat <- object$modMat) # number of _multivariate_
                                          # observations
    nAxes <- addArgs$nAxes # number of axes or (latent dimensions)
    trmDims <- c(nVar = nVar, nObs = nObs, nAxes = nAxes)

    if(is.null(addArgs$seed)) {
        obsInds <- apply(obsMat == 1, 1, which)
        obsFac <- as.factor(colnames(obsMat)[obsInds])
        initLoad <- svd(2 * xtabs(auxEnv$response ~ obsFac + varFac) - 1)$v[,1:nAxes]
        loadMat <- repSparseFull(nVar, nAxes, as.vector(initLoad))
    } else {
        set.seed(addArgs$seed)
        loadMat <- repSparseFull(nVar, nAxes, rnorm(nVar * nAxes))
    }
    
    initLoadMat <- as.matrix(loadMat)
    ut <- upper.tri(initLoadMat, diag = FALSE)
    initLoadMat[ut] <- 0
    loadMat$trans <- local({
        trmDims <- trmDims
        initLoadMat <- initLoadMat
        lt <- lower.tri(initLoadMat, diag = TRUE)
        init <- initLoadMat[lt]
        function(matPars) {
            initLoadMat[lt] <- matPars
            as.numeric(initLoadMat)
        }
    })
    loadMat <- update(loadMat)
    
    Zt <- kr(t(subset(loadMat, as.numeric(varFac))),
             t(resetTransConst(simplifyRepSparse(as.repSparse(obsMat)))))

    ## check equivalence of subset versus matrix multiplication
    ## plot(t(t(as.matrix(loadMat)) %*% as.matrix(as.repSparse(obsMat))),
    ##      as.matrix(obsMat))

    Lambdat <- resetTransConst(repSparseIdent(nrow(Zt)))
    
    lowerLoads <- rep(-Inf, length(getInit(Zt)))
    upperLoads <- rep( Inf, length(getInit(Zt)))

    packReTrm(object, Zt, Lambdat,
              lowerLoads = lowerLoads,
              upperLoads = upperLoads,
              devfunEnv = devfunEnv)
}

##' Random effects term for structural equation model
##'
##' @export
##' @template setReTrm
##' @templateVar cls sem
##' @templateVar form \code{sem(1 | grpFac, loadMat)}
##' @templateVar arg c("grpFac", "nAxes", "seed")
##' @templateVar desc c("grouping factor for loadings", "number of axes", "random seed for initial axis values")
setReTrm.sem <- function(object, addArgsList,
                         auxEnv = NULL, devfunEnv = NULL) {
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    nl <- nlevels(grpFac <- object$grpFac)
    nc <- addArgs$nAxes
    set.seed(addArgs$seed)
    ## loadMat <- rRepSparse(nl, nc, nl * nc, nl * nc)
    modMat <- subset(loadMat, as.numeric(grpFac))
    indMat <- resetTransConst(as.repSparse(addArgs$obsFac))

    ## check equivalence of subset versus matrix multiplication
    ## plot(t(t(as.matrix(loadMat)) %*% as.matrix(as.repSparse(grpFac))),
    ##      as.matrix(modMat))

    Zt <- kr(t(modMat), indMat) ## FIXME: obsFac before modMod??
    Lambdat <- resetTransConst(repSparseIdent(nrow(Zt)))
    
    lowerLoads <- rep(-Inf, length(getInit(Zt)))
    upperLoads <- rep( Inf, length(getInit(Zt)))

    packReTrm(object, Zt, Lambdat,
              lowerLoads = lowerLoads,
              upperLoads = upperLoads,
              devfunEnv = devfunEnv)
}



##' Random effects term for covariance proportional to identity matrix
##'
##' @export
##' @template setReTrm
##' @templateVar cls identity
##' @templateVar form \code{identity(linForm | grpFac)}
##' @templateVar arg c("linForm", "grpFac")
##' @templateVar desc c("linear model formula decribing effects", "grouping factor")
setReTrm.identity <- function(object, addArgsList,
                              auxEnv = NULL, devfunEnv = NULL) {

                                        # Zt
    Zt <- resetTransConst(kr(as.repSparse(object$grpFac),
                           t(as.repSparse(object$modMat))))

                                        # Lambdat
    nl <- nlevels(object$grpFac)
    nc <-    ncol(object$modMat)
    Lambdat <- repSparseIdent(nl * nc)

                                        # pack
    packReTrm(object, Zt, Lambdat, devfunEnv = devfunEnv)
}

##' Random effects term for a flexible variance function
##'
##' @export
##' @template setReTrm
##' @templateVar cls flexvar
##' @templateVar form \code{flexvar(linForm, init, nBasis)}
##' @templateVar arg c("linForm", "init", "nBasis")
##' @templateVar desc c("linear model formula decribing effects", "initial values for coefficients for the linear predictor of the variance function", "number of coefficients for the linear predictor")
setReTrm.flexvar <- function(object, addArgsList,
                             auxEnv = NULL, devfunEnv = NULL) {
                                        # get additional arguments and
                                        # sample size
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    n <- nrow(object$modMat)

                                        # Zt
    Zt <- resetTransConst(repSparseIdent(n))

                                        # Lambdat
    inds <- seq_len(n); baselineVars <- rep(1, n)
    if(is.null(init <- addArgs$init)) init <- rep(0, addArgs$nBasis)
    Lambdat       <- repSparseDiag  (      baselineVars, inds     )
    Lambdat$trans <- mkFlexDiagTrans(init, baselineVars, devfunEnv)

                                        # pack
    packReTrm(object, Zt, Lambdat,
              lowerCovar = rep(-Inf, length(init)),
              devfunEnv = devfunEnv)
}


##' Random effects term with exponential decay in distance-based covariance
##'
##' @export
##' @template setReTrm
##' @templateVar cls expDecay
##' @templateVar form \code{expDecay(1 | grpFac, distCutoff, minCov, distMat)}
##' @templateVar arg c("grpFac", "distCutoff", "minCov", "distMat")
##' @templateVar desc c("grouping factor (e.g. with levels given by geographical sites)", "maximum distance with covariance greater than minCov", "minimum covariance", "distance matrix object over the levels of grpFac")
setReTrm.expDecay <- function(object, addArgsList,
                              auxEnv = NULL, devfunEnv = NULL) {
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    distCutoff         <- addArgs$distCutoff
    minCov             <- addArgs$minCov
    distMat <- edgeMat <- addArgs$distMat
    init               <- addArgs$init

    if(is.null(init)) init <- c(1, 1)

                                        # approximate exponential
                                        # decay model by an edge-based
                                        # model, with edges between
                                        # levels that are less than
                                        # distCutoff from each other
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
    Zt <- resetTransConst(kr(Jt, t(as.repSparse(object$modMat))))
    
    transFun <- local({
        init <- init
        edgeDists <- edgeDists
        distCutoff <- distCutoff
        minCov <- minCov
        function(matPars) {
            q1 <- (minCov - 1)/(exp(-(matPars[1]) * distCutoff) - 1)
            q2 <- 1 - q1
            matPars[2] * (q2 + q1 * exp(-(matPars[1]) * edgeDists))
        }
    })
    Lambdat <- repSparseDiag(transFun(1))
    Lambdat$trans <- transFun
    packReTrm(object, Zt, update(Lambdat), devfunEnv = devfunEnv)
}

##' Random effects term for edge-based model of tree-induced covariance
##'
##' @export
##' @template setReTrm
##' @templateVar cls phyloEdge
##' @templateVar form \code{edge(linForm | grpFac, phy)}
##' @templateVar arg c("linForm", "grpFac", "phy")
##' @templateVar desc c("linear model formula decribing effects", "grouping factor", "phylo object relating the levels of the grouping factor")
setReTrm.phyloEdge <- function(object, addArgsList,
                               auxEnv = NULL, devfunEnv = NULL) {
                                        # get additional arguments
    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)

                                        # Zt
    Jedge <- as(edgeTipIndicator(addArgs$phy), "sparseMatrix")
    Jspp <- as(object$grpFac, "sparseMatrix")
    Jt <- resetTransConst(as.repSparse(Jedge %*% Jspp))
    Zt <- resetTransConst(kr(Jt, t(as.repSparse(object$modMat))))

                                        # Lambdat
    nl <- nrow(Jedge)
    nc <- ncol(object$modMat)
    Lambdat <- repSparseIdent(nl * nc)
    
                                        # pack
    packReTrm(object, Zt, Lambdat,
              devfunEnv = devfunEnv)
}

##' Random effects term for cooccurence model (experimental)
##'
##' @export
##' @template setReTrm
##' @templateVar cls cooccur
##' @templateVar form \code{cooccur(linForm | grpFac)}
##' @templateVar arg c("linForm", "grpFac")
##' @templateVar desc c("linear model formula decribing effects", "grouping factor")
setReTrm.cooccur <- function(object, addArgsList,
                             auxEnv = NULL, devfunEnv = NULL) {

                                        # Zt
    Jt <- as.repSparse(as(object$grpFac, "sparseMatrix"))
    Zt <- resetTransConst(kr(Jt, t(as.repSparse(object$modMat))))

                                        # Lambdat
    nCovPars <- choose(ncol(object$modMat), 2)
    Tt <- t(repSparseCorMatChol(rep(0, nCovPars)))
    Lambdat <- rep(Tt, nrow(Jt), type = "diag")
    
                                        # pack
    packReTrm(object, Zt, Lambdat,
              devfunEnv = devfunEnv)
}

##' Random effects term with constant variance within groups
##'
##' @export
##' @template setReTrm
##' @templateVar cls varIdent
##' @templateVar form \code{varIdent(1 | grpFac)}
##' @templateVar arg c("grpFac")
##' @templateVar desc c("grouping factor")
setReTrm.varIdent <- function(object, addArgsList,
                              auxEnv = NULL, devfunEnv = NULL) {
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
              lowerCovar = rep(0, nl),
              devfunEnv = devfunEnv)
}

##' Random effects term with exponential variance structure
##'
##' @export
##' @template setReTrm
##' @templateVar cls varExp
##' @templateVar form \code{varExp(linForm | grpFac)} or \code{varExp(linForm)}
##' @templateVar arg c("linForm", "grpFac")
##' @templateVar desc c("linear model formula decribing effects", "grouping factor")
setReTrm.varExp <- function(object, addArgsList,
                            auxEnv = NULL, devfunEnv = NULL) {
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
              lowerCovar = rep(-Inf, length(init)),
              devfunEnv = devfunEnv)
}

##' Random effects term with observation level random effects
##'
##' @export
##' @template setReTrm
##' @templateVar cls obslev
##' @templateVar form \code{obslev(1)}
##' @templateVar arg "no arguments"
##' @templateVar desc "NULL"
setReTrm.obslev <- function(object, addArgsList,
                             auxEnv = NULL, devfunEnv = NULL) {

    n <- nrow(object$modMat)
    Zt <- resetTransConst(repSparseIdent(n))
    Lambdat <- repSparseIdent(n)
    packReTrm(object, Zt, Lambdat,
              devfunEnv = devfunEnv)
}

##' Random effects term from an \code{nlme}-style \code{corStruct}
##' object
##'
##' @export
##' @template setReTrm
##' @templateVar cls nlmeCorStruct
##' @templateVar form \code{nlmeCorStruct(1, corObj)} or \code{nlmeCorStruct(1 | grpFac, corObj)}
##' @templateVar arg c("corObj", "sig")
##' @templateVar desc c("corStruct object", "initial standard deviation")
setReTrm.nlmeCorStruct <- function(object, addArgsList,
                                   auxEnv = NULL, devfunEnv = NULL) {

    addArgs <- getAddArgs(object$addArgs[-1], addArgsList)
    corObj <- addArgs$corObj
    sig <- addArgs$sig

    modMat <- object$modMat
    if(is.na(object$grpFac[[1]])) {
        grpFac <- as.factor(1:nrow(modMat))
    } else {
        grpFac <- object$grpFac
    }
    Lambdat <- repSparseCorFactor(corObj, sig = sig)
    Zt <- resetTransConst(kr(as.repSparse(grpFac),
                           t(as.repSparse(modMat))))
    packReTrm(object, Zt, Lambdat,
              devfunEnv = devfunEnv)
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
                      lowerCovar, upperCovar,
                      devfunEnv = NULL) {
    if(missing(lowerLoads)) lowerLoads <- setLowerDefault(getInit(Zt))
    if(missing(upperLoads)) upperLoads <- setUpperDefault(getInit(Zt))
    if(missing(lowerCovar)) lowerCovar <- setLowerDefault(getInit(Lambdat))
    if(missing(upperCovar)) upperCovar <- setUpperDefault(getInit(Lambdat))
    structure(c(object,
                lme4:::namedList(Zt, Lambdat,
                                 lowerLoads, upperLoads,
                                 lowerCovar, upperCovar,
                                 devfunEnv)),
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

##' @export
VarCorr.factAnal <- function(x, sigma = 1, rdig = 3) {
    tcrossprod(loadings.factAnal(x))
}

##' Update a random effects term structure with new parameters
##'
##' @param object a \code{reTrmStruct} object
##' @param newCovar new covariance parameters
##' @param newLoads new loadings parameters
##' @param ... potential additional arguments (ignored by the default
##' method)
##' @rdname update.reTrmStruct
##' @method update reTrmStruct
##' @export
update.reTrmStruct <- function(object, newCovar, newLoads, ...) {
    object$Lambdat <- update(object$Lambdat, newCovar)
    object$Zt      <- update(object$Zt,      newLoads)
    return(object)
}

##' @rdname update.reTrmStruct
##' @method update flexvar
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


## update.factAnal <- function(object, newCovar, newLoads, ...) {
##     object <- update.reTrmStruct(object, newCovar, newLoads)
##     transEnv <- environment(object$Lambdat$trans)
##     assignWith(expr  = resp$y,
##                name  = "respVar",
##                data  = transEnv$devfunEnv,
##                envir = transEnv)
##     return(object)
## }

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
simAddArgs.phyloEdge <- function(object, rtreeArgs = list(),
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

##' Simulate from a random effects term
##'
##' @param object \code{reTrmStruct} object
##' @export
simReTrm <- function(object) {
    with(object, {
        as.numeric(rnorm(nrow(Lambdat)) %*% as.matrix(Lambdat) %*% as.matrix(Zt))
    })
}

##' @param type see \code{\link{ranef.strucGlmer}}
##' @param ... not used (for consistency with generic)
##' @rdname setReTrm
##' @export
ranef.reTrmStruct <- function(object, type = c("u", "Lu", "ZLu"), ...) {
    type <- type[1]
    trmName <- paste(object$grpName, class(object)[[1]], sep = ".")
    pp <- object$devfunEnv$pp
    nms <- names(nRePerTrm <- object$devfunEnv$nRePerTrm)
    if (type == "u") {
        re <- pp$u(1)
    } else {
        re <- pp$b(1)
        if (type == "ZLu") b <- subRagByLens(re, nRePerTrm)[[trmName]]
        return(as.matrix(t(object$Zt), sparse = TRUE %*% b))
    }
    subRagByLens(re, nRePerTrm)[[which(nms == trmName)]]
}

##' @param x \code{\link{setReTrm.factAnal}} object
##' @param ... not used (for consistency with generic in \code{vegan}
##' package)
##' @importFrom vegan scores
##' @rdname setReTrm.factAnal
##' @export
scores.factAnal <- function(x, ...) {
    ## requires vegan
    trans <- environment(x$Zt$trans)$Btrans
    trmDims <- environment(trans)$trmDims
    trmLoads <- matrix(trans(getInit(x$Zt)),
                       trmDims["nVar"],
                       trmDims["nAxes"])
    trmFacts <- matrix(ranef(x),
                       trmDims["nObs"],
                       trmDims["nAxes"])
    colnames(trmLoads) <- colnames(trmFacts) <- paste("axis", 1:trmDims["nAxes"])
    list(loadings = trmLoads, factors = trmFacts)
}

##' @importFrom stats biplot
##' @rdname setReTrm.factAnal
##' @export
biplot.factAnal <- function(x, ...) {
    l... <- list(...)
    do.call(biplot, c(setNames(scores(x), c("y", "x")), l...))
}

##' Bind repeated sparse matrices and sort their indices to be
##' compatible with column-compressed sparse matrices
##'
##' @param mats list of \code{repSparse} matrix objects
##' @param type type of bind
sortedBind <- function(mats, type = c("row", "col", "diag")) {
    standardSort(.bind(mats, type = type))
}

##' Make transformation function to be used with a column-compressed
##' sparse matrix
##'
##' @param object repeated sparse matrix object
##' @export
##' @examples
##' set.seed(1)
##' X <- rRepSparse(3, 7, 2, 12)
##' (Xsparse <- as.matrix(X, sparse = TRUE))
##' Xtrans <- mkSparseTrans(X)
##' slot(Xsparse, "x") <- Xtrans(rnorm(2))
##' print(Xsparse)
mkSparseTrans <- function(object) {
    ans <- function(matPars) trans(matPars)[inds]
    environment(ans) <- new.env(parent = environment(object$trans))
    environment(ans)$trans <- object$trans
    environment(ans)$inds <- object$valInds
    return(ans)
}

