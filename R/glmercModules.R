##' Parse mixed model formula for levels-covariance model
##'
##' @param formula mixed model formula
##' @param data data
##' @param family family
##' @param covList list of covariance matrices for random effects
##' grouping factors in \code{formula}
##' @param strList list of structure matrices for random effects
##' grouping factors in \code{formula}
##' @param tileCov in kronecker products, should the covariance
##' matrices in \code{covList} be tiled (\code{tileCov = TRUE}) or
##' distributed (\code{tileCov = FALSE})?
##' @param giveCsparse use compressed-form \code{CsparseMatrix}
##' representations (\code{TRUE}), or triplet-form
##' \code{TsparseMatrix} representations (\code{FALSE})?
##' @param ... not used yet
##' @importFrom lme4 findbars nobars
##' @export
glmercFormula <- function(formula, data = NULL, family = binomial,
                          covList, strList,
                          tileCov = TRUE,
                          giveCsparse = TRUE, ...) {

                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

                                        # append the covariance and
                                        # structure matrices over the
                                        # levels associated with each
                                        # random effects term
    for(i in seq_along(reTrmsList)) {
        reTrmsList[[i]]$strMat <- strList[[reTrmsList[[i]]$grpName]]
        if(is.null(reTrmsList[[i]]$strMat)) { # if no str mat, use
                                              # identify
            nli <- nlevels(reTrmsList[[i]]$grpFac)
            reTrmsList[[i]]$strMat <- diag(1, nli, nli)
        }
        reTrmsList[[i]]$covMat <- covList[[reTrmsList[[i]]$grpName]]
        if(is.null(reTrmsList[[i]]$covMat)) { # if no cov mat, use
                                              # identity
            nstr <- nrow(reTrmsList[[i]]$strMat)
            reTrmsList[[i]]$covMat <- diag(1, nstr, nstr)
        }
    }

                                        # extract the model matrices,
                                        # grouping factors, and cov
                                        # mats
    modMats <- lapply(reTrmsList, "[[", "modMat")
    grpFacs <- lapply(reTrmsList, "[[", "grpFac")
    covMats <- lapply(reTrmsList, "[[", "covMat")
    strMats <- lapply(reTrmsList, "[[", "strMat")

                                        # set up constant grouping
                                        # factors and 1 by 1 cov mats
                                        # for mkTemplateReTrms (this
                                        # is where we could add
                                        # nesting)
    grpFacConst <- rep(list(as.factor(rep("const", nrow(data)))), length(grpFacs))
    covMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))
    strMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))

                                        # construct the full random
                                        # effects model matrix,
                                        # relative covariance factor,
                                        # etc.
    if(tileCov) {
        re <- mkTemplateReTrms(modMats,
                               grpFacs, grpFacConst,
                               covMats, covMatConst,
                               strMats, strMatConst,
                               giveCsparse)
    } else {
        re <- mkTemplateReTrms(modMats,
                               grpFacConst, grpFacs,
                               covMatConst, covMats,
                               strMatConst, strMats,
                               giveCsparse)
    }

    
                                        # organize output value
    return(c(re,
             list(X = model.matrix(nobars(formula), data),
                  y = model.response(model.frame(nobars(formula), data)),
                  flist = simplifyFacList(grpFacs),
                  modMats = modMats,
                  strMats = strMats)))
}

##' Get components of a \code{glmerc} object
##'
##' @param object \code{\link{glmerc}} object
##' @param name object name
##' @export
getMEc <- function(object, name = c("X", "Z", "Zt", "y", "TmodMat", "cnms")) {
    if (missing(name)) 
        stop("'name' must not be missing")
    if (length(name <- as.character(name)) > 1) {
        names(name) <- name
        return(lapply(name, getMEc, object = object))
    }
    cnms <- lapply(object$parsedForm$modMats, colnames)
    TmodMat <- object$parsedForm$TmodMat
    
    for(i in 1:length(TmodMat)) {
        TmodMat[[i]]@x <- covarByTerms(object)[[i]]
        dimnames(TmodMat[[i]]) <- rep(cnms[i], 2)
    }
    switch(name,
           X = object$X,
           Z = t(object$parsedForm$Zt),
           Zt = object$parsedForm$Zt,
           y = object$y,
           TmodMat = TmodMat,
           cnms = cnms)
}

##' Make covariance template random effects term
##'
##' @param modMat a model matrix
##' @param grpFac1,grpFac2 grouping factors (typically dimension IDs
##' of a melted response matrix)
##' @param covMat1,covMat2 covariance matrices among the levels of
##' \code{grpFac1} and \code{grpFac2}
##' @param strMat1,strMat2 structure matrices
##' @param giveCsparse use compressed-form \code{CsparseMatrix}
##' representations (\code{TRUE}), or triplet-form
##' \code{TsparseMatrix} representations (\code{FALSE})?
##' @export
mkTemplateReTrm <- function(modMat,
                            grpFac1, grpFac2,
                            covMat1, covMat2,
                            strMat1, strMat2,
                            giveCsparse = TRUE) {

    badDims <- FALSE # FIXME: improve condition
        # (length(levels(grpFac1)) != nrow(covMat1)) ||
        # (length(levels(grpFac2)) != nrow(covMat2))
    if(badDims) stop("covariance matrix incompatible with its grouping factor")

                                        # use consistent Matrix
                                        # classes (Triplet-form
                                        # TSparseMatrix -- i, j, x --
                                        # not Compressed-form
                                        # CSparseMatrix)
    matClass <- if(giveCsparse) "CsparseMatrix" else "TsparseMatrix"

                                        # expand model matrix from
                                        # indicator matrices (double
                                        # `as` required because
                                        # currently no method in
                                        # Matrix to go from factor to
                                        # TsparseMatrix) (FIXME: make
                                        # sure levels of grpFac's are
                                        # in the same order as the
                                        # dimnames of the covMat's)
    J1 <- as(strMat1 %*% as(grpFac1, "sparseMatrix"), matClass)
    J2 <- as(strMat2 %*% as(grpFac2, "sparseMatrix"), matClass)
    Zt <- KhatriRao(KhatriRao(J2, t(modMat)), J1)

    nc <- ncol(modMat)
    rowIndices <- sequence(1:nc) # note different order from lme4
                                 # (FIXME: necessary?)
    colIndices <- rep(1:nc, 1:nc)
    TmodMat <- TmodMatLind <- TmodMatOnes <-
        sparseMatrix(rowIndices, colIndices,
                     x = 1 * (rowIndices == colIndices),
                     giveCsparse = giveCsparse)
    TmodMatLind@x <- as.numeric(seq_along(TmodMat@x))
    TmodMatOnes@x <- as.numeric(rep(1, length(TmodMat@x)))

    T1 <- T1ones <- as(chol(covMat1), matClass)
    T2 <- T2ones <- as(chol(covMat2), matClass)
    T1ones@x <- as.numeric(rep(1, length(T1ones@x)))
    T2ones@x <- as.numeric(rep(1, length(T2ones@x)))

    LambdatBaseline <- T2 %x% TmodMatOnes %x% T1
    LambdatLind <- T2ones %x% TmodMatLind %x% T1ones
    LambdatCovar <- T2ones %x% TmodMat %x% T1ones
    Lambdat <- T2 %x% TmodMat %x% T1

    return(list(Zt = Zt,
                Lambdat = Lambdat, 
                LambdatLind = LambdatLind,
                LambdatCovar = LambdatCovar,
                LambdatBaseline = LambdatBaseline,
                TmodMat = TmodMat,
                nCovar = length(TmodMat@x)))
}

##' Transpose a list
##'
##' @param lst \code{\link{list}}
##' @export
listTranspose <- function(lst) {
    lstExtract <- function(i) lapply(lst, "[[", i)
    setNames(lapply(names(lst[[1]]), lstExtract), names(lst[[1]]))
}

reTrmsBdiag <- function(lst) {
    ii <- unlist(lapply(lst, slot, "i"))
    jj <- unlist(lapply(lst, slot, "j"))
    xx <- unlist(lapply(lst, slot, "x"))
    sparseMatrix(i = ii, j = jj, x = xx, giveCsparse = FALSE)
}


##' Make several covariance template random effects terms
##'
##' @param modMat,grpFac1,grpFac2,covMat1,covMat2,strMat1,strMat2,giveCsparse lists of inputs to
##' \code{\link{mkTemplateReTrm}}
##' @export
mkTemplateReTrms <- function(modMat,
                             grpFac1, grpFac2,
                             covMat1, covMat2,
                             strMat1, strMat2,
                             giveCsparse) {

    reTrmsList <- listTranspose(mapply(mkTemplateReTrm,
                                       modMat,
                                       grpFac1, grpFac2,
                                       covMat1, covMat2,
                                       strMat1, strMat2,
                                       SIMPLIFY = FALSE,
                                       MoreArgs = list(giveCsparse = giveCsparse)))

    if(!(length(reTrmsList$Zt) == 1L)) {
        for(i in 2:length(reTrmsList$LambdatLind)) {
            j <- i - 1
            reTrmsList$LambdatLind[[i]]@x <-
                reTrmsList$LambdatLind[[i]]@x + max(reTrmsList$LambdatLind[[j]]@x)
        }
    }
    within(reTrmsList, {
        Zt <- do.call(rBind, Zt)
        
        Lambdat <- .bdiag(Lambdat)
        LambdatLind <- .bdiag(LambdatLind)
        LambdatCovar <- .bdiag(LambdatCovar)
        LambdatBaseline <- .bdiag(LambdatBaseline)
        
        Lind <- LambdatLind@x
        covar <- LambdatCovar@x[match(sort(unique(Lind)), Lind)]
        baseline <- LambdatBaseline@x
        mapToCovFact <- function(covar) baseline * covar[Lind]
    })
}

##' Make covariance factor for a covariance template term
##'
##' @param covMat template covariance matrix
##' @param nLevels number of levels over which the grouping factor varies
##' @param diagModel templates on the block diagonal (\code{TRUE}) or
##' densely tiled over the entire matrix (\code{FALSE})?
mkTemplateTermLambdat <- function(covMat, nLevels, diagModel = TRUE) {
    ## cm <- crossprod(matrix(rnorm(25), 5, 5))
    ## LamtDiag <- mkTemplateTermLambdat(cm, 4, TRUE)
    ## LamtDense <- mkTemplateTermLambdat(cm, 4, FALSE)
    ## image(LamtDiag)
    ## image(LamtDense)
    Lt <- as(chol(covMat), "sparseMatrix")
    if(diagModel) {
        return(.bdiag(rep(list(Lt), nLevels)))
    } else {
        n <- nrow(covMat)
        zeros <- matrix(0, (nLevels - 1) * n, nLevels * n)
        Lts <- do.call(cBind, rep(list(Lt), nLevels))
        return(rBind(Lts, zeros))
    }
}

mkTemplateTermZt <- function(explVar, grpFac) {
    J <- as(as.factor(grp), "sparseMatrix")
    explVar 
}

##' Make function for computing deviance profiles
##' 
##' @param mod \code{\link{glmerc}} object
##' @param whichPar which parameter to profile
##' @return function for computing the profiled deviance over a single
##' parameter
##' @export
mkPfun <- function(mod, whichPar) {
    local({
        dfun <- mkNfun(mod)
        optPar <- mod$opt$par
        whichPar <- whichPar
        op <- optPar[-whichPar]
        function(par) {
            fn <- function(op) {
                parNow <- numeric(length(op) + 1)
                parNow[whichPar] <- par
                parNow[-whichPar] <- op
                dfun(parNow)
            }
            opt <- minqa:::bobyqa(op, fn, lower = mod$lower[-whichPar],
                          control = list(iprint = 0L, maxfun = 500))
            op <<- opt$par # FIXME: assign to local env
            return(opt$fval)
        }
    })
}

##' Make function for computing deviance slices
##' 
##' @param mod \code{\link{glmerc}} object
##' @param whichPar which parameter to slice
##' @return function for computing 1D slices through the deviance at
##' the optimum
##' @export
mkSfun <- function(mod, whichPar) {
    local({
        dfun <- mkNfun(mod)
        whichPar <- whichPar
        par <- mod$opt$par
        function(par1D) {
            par[whichPar] <- par1D
            return(dfun(par))
        }
    })
}

##' Make function for computing normalized deviance
##'
##' @param mod model object
##' @return function for computing the normalized deviance, which
##' equals zero at the optimum
##' @export
mkNfun <- function(mod) {
    local({
        dfun <- mod$dfun
        optFunVal <- dfun(mod$opt$par)
        function(par) dfun(par) - optFunVal
    })
}
