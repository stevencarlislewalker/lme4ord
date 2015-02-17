##' Make covariance template random effects term
##'
##' @param modMat a model matrix
##' @param grpFac1,grpFac2 grouping factors (typically dimension IDs
##' of a melted response matrix)
##' @param covMat1,covMat2 covariance matrices among the levels of
##' \code{grpFac1} and \code{grpFac2}
##' @export
mkTemplateReTrm <- function(modMat, grpFac1, grpFac2, covMat1, covMat2) {

    badDims <-
        (length(levels(grpFac1)) != nrow(covMat1)) ||
        (length(levels(grpFac2)) != nrow(covMat2))
    if(badDims) stop("covariance matrix incompatible with its grouping factor")

                                        # use consistent Matrix
                                        # classes (Triplet-form Sparse
                                        # Matrix -- i, j, x -- not
                                        # Compressed-form Sparse
                                        # Matrix)
    matClass <- "TsparseMatrix"
    giveCsparse <- FALSE

                                        # expand model matrix from
                                        # indicator matrices (double
                                        # `as` required because
                                        # currently no method in
                                        # Matrix to go from factor to
                                        # TsparseMatrix) (FIXME: make
                                        # sure levels of grpFac's are
                                        # in the same order as the
                                        # dimnames of the covMat's)
    J1 <- as(as(grpFac1, "sparseMatrix"), matClass)
    J2 <- as(as(grpFac2, "sparseMatrix"), matClass)
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
    LambdatTheta <- T2ones %x% TmodMat %x% T1ones
    Lambdat <- T2 %x% TmodMat %x% T1
    
    return(list(Zt = Zt,
                Lambdat = Lambdat, 
                LambdatLind = LambdatLind,
                LambdatTheta = LambdatTheta,
                LambdatBaseline = LambdatBaseline))
}


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
##' @param modMat,grpFac1,grpFac2,covMat1,covMat2 lists of inputs to
##' \code{\link{mkTemplateReTrm}}
##' @export
mkTemplateReTrms <- function(modMat, grpFac1, grpFac2, covMat1, covMat2) {

    reTrmsList <- listTranspose(mapply(mkTemplateReTrm,
                                       modMat, grpFac1, grpFac2, covMat1, covMat2,
                                       SIMPLIFY = FALSE))

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
        LambdatTheta <- .bdiag(LambdatTheta)
        LambdatBaseline <- .bdiag(LambdatBaseline)
        
        Lind <- LambdatLind@x
        theta <- LambdatTheta@x[match(sort(unique(Lind)), Lind)]
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
##' @examples
##' cm <- crossprod(matrix(rnorm(25), 5, 5))
##' LamtDiag <- mkTemplateTermLambdat(cm, 4, TRUE)
##' LamtDense <- mkTemplateTermLambdat(cm, 4, FALSE)
##' image(LamtDiag)
##' image(LamtDense)
mkTemplateTermLambdat <- function(covMat, nLevels, diagModel = TRUE) {
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

