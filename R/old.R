## # ----------------------------------------------------------------------
## # factorCov
## # ----------------------------------------------------------------------

## ## Factors with covariance matrices over the levels
## ##
## ## @param fac a factor
## ## @param covMat a covariance matrix over the levels (if missing use
## ## identity matrix)
## ## @param stan standardize covariance matrix to determinant one?
## ## @return a \code{factorCov} object
## ## @rdname factorCov
## ## @export
## factorCov <- function(fac, covMat, stan = TRUE) {
##     fac <- as.factor(fac)
##     if(missing(covMat)) {
##         covMat <- diag(1, nlevels(fac), nlevels(fac))
##     } else if(stan) {
##         covMat <- stanCov(covMat)
##     }
##     if(!(nlevels(fac) == nrow(covMat))) stop("covMat must be over levels(fac)")
##     attr(fac, "covMat") <- covMat
##     class(fac) <- c("factorCov", "factor")
##     return(fac)
## }

## ## @rdname factorCov
## ## @export
## print.factorCov <- function(x, ...) {
##     ret <- x
##     attr(x, "covMat") <- NULL
##     class(x) <- "factor"
##     print(as.factor(x), ...)
##     cat("Note: this factor includes a covariance matrix over its levels, attr(., \"covMat\")\n")
## }


## # ----------------------------------------------------------------------
## # logisticPCA
## # ----------------------------------------------------------------------

## ## Construct deviance function for logistic principal components analysis
## ##
## ## @param Y response matrix
## ## @param d number of axes
## ## @param Xrow row model matrix
## ## @param Xcol col model matrix
## ## @param Upenalty penalty for axis scores
## ## @param thetaPenalty penalty for theta for the axes
## ## @param Ustart optional initial matrix of column scores
## ## @param thetaStart optional initial value for theta vector
## ## @importFrom Matrix t Diagonal rBind sparseMatrix .bdiag bdiag
## ## @export
## logisticPcaDevfun <-
## function(Y, d, Xrow, Xcol, Ustart, thetaStart, Upenalty = 0, thetaPenalty = 0) {

##     y <- as.numeric(Y)

##     t <- ncol(Y)
##     s <- nrow(Y)
##     n <- t*s
    
##     if(missing(Ustart)) U <- initU(Y, d) else U <- Ustart
##     Uind <- U

##     if(missing(Xrow)) 1
##     if(missing(Xcol)) 1
        
    
##     ## subscript explanations:
##     ## 
##     ## U -- based on U matrix
##     ## s -- related to sites
##     ## t -- related to types (e.g. species)
##     ## O -- based on Other possible random effects terms
##     ZtAxes <- kronecker(t(U), Diagonal(s))
##     JtRow <- as(factor(rep(1:s, times = t), 1:s), "sparseMatrix")
##     JtCol <- as(factor(rep(1:t,  each = s), 1:t), "sparseMatrix")
##     ## Xt_s <- rbind(1, rep(x, times = t))
##     ## Xt_t <- matrix(1, ncol = n)
##     ## Xt_t <- rbind(1, rep(z,  each = s))
##     ## ZtIntrcp <- rBind(KhatriRao(JtRow, Xt_s),
##     ##                   KhatriRao(JtCol, Xt_t))
##     ## Xt_t <- matrix(1, ncol = n)
##     ## Xt_s <- matrix(1, ncol = n)
##     XIntrcp <- matrix(1, ncol = n)
##     Zt_O <- rBind(JtRow, JtCol)
##     X <- t(XIntrcp)
    
##     Zt <- rBind(ZtAxes, Zt_O)

##     Zmap <- local({
##                                         # point to 
##         Uind@x <- as.numeric(1:length(Uind@x))
##         Zind_U <- kronecker(t(Uind), Diagonal(s))
##         Zind_O <- Zt_O
##         Zind_O@x <- rep(Inf, length(Zind_O@x))
##         Zind <- rBind(Zind_U, Zind_O)@x
##         whichToUpdate <- is.finite(Zind)
##         ZindUpdate <- Zind[whichToUpdate]
##         out <- Zt@x
##         function(phi) {
##             out[whichToUpdate] <- phi[ZindUpdate]
##             return(as.numeric(out))
##         }
##     })
    

##     #phi <- round(rnorm(length(U@x)), 2)
##     #Zt@x <- Zmap(phi)
##     phi <- U@x
##     Zt@x <- Zmap(phi)

##     Lambdat_U <- Diagonal(d*s)
##     theta_U <- 1 # rep(1, d)

##     nc_s <- 1 # nrow(Xt_s)
##     tempInd_s <- 1:nc_s
##     ii_s <- rep(tempInd_s, tempInd_s)
##     jj_s <- sequence(tempInd_s)
##     temp_s <- sparseMatrix(i = ii_s, j = jj_s,
##                            x = 1*(ii_s==jj_s))
##     theta_s <- temp_s@x
##     Lambdat_s <- .bdiag(rep(list(t(temp_s)), s))

##     nc_t <- 1 # nrow(Xt_t)
##     tempInd_t <- 1:nc_t
##     ii_t <- rep(tempInd_t, tempInd_t)
##     jj_t <- sequence(tempInd_t)
##     temp_t <- sparseMatrix(i = ii_t, j = jj_t,
##                            x = 1*(ii_t==jj_t))
##     theta_t <- temp_t@x
##     Lambdat_t <- .bdiag(rep(list(t(temp_t)), t))
##     Lambdat_O <- bdiag(Lambdat_s, Lambdat_t)
##     Lambdat <- bdiag(Lambdat_U, Lambdat_O)

##     Lind_U <- rep(1, d*s) # rep(1:d, each = s)
##     Lind_s <- rep(ii_s + nc_s*(jj_s - 1) - choose(jj_s, 2), s)
##     Lind_t <- rep(ii_t + nc_t*(jj_t - 1) - choose(jj_t, 2), t)
##     Lind_O <- c(Lind_s, Lind_t + max(Lind_s))
##     Lind <- c(Lind_U, Lind_O + max(Lind_U))

##     theta_O <- c(theta_s, theta_t)
##     theta <- c(theta_U, theta_O)
##     if(!missing(thetaStart)) theta[] <- thetaStart
    
##     resp <- new(Class = "glmResp",
##                 family = binomial(link="logit"),
##                 y = y)
##                                         # test if on flexLambda branch
##                                         # of lme4.  TODO: eliminate
##                                         # and warn once flexLambda is
##                                         # on CRAN
##     ifl <- !is.character(try(lme4::isFlexLambda(), TRUE)) 
##     if(ifl){
##         pp <- new(Class = "merPredD",
##                   X = X,
##                   Zt = Zt,
##                   Lambdat = Lambdat,
##                   Lind = Lind,
##                   thfun = local({Lind <- Lind; function(theta) theta[Lind]}),
##                   theta = theta,
##                   phi = numeric(0),
##                   phifun1 = function(phi) { }, ## no-op
##                   n = n)
##     } else {
##         pp <- new(Class = "merPredD",
##                   X = X,
##                   Zt = Zt,
##                   Lambdat = Lambdat,
##                   Lind = Lind,
##                   theta = theta,
##                   n = n)
##     }
##     function(pars) {    
##         resp$updateMu(pp$linPred(1) + 0)
##                                         # zth are the elements of U
##         zth <- pars[-(1:length(theta[-1]))]
##         zx <- Zmap(zth)
##                                         # th are the elements of Lambda
##         th <- c(1/mean(abs(zth)), pars[1:length(theta[-1])])

##         pp$setTheta(as.double(th))
##         pp$setZt(as.double(zx))
##         p <- glmerLaplaceHandle(pp$ptr(),
##                                 resp$ptr(),
##                                 0, 1e-3, 30, 0)
##         resp$updateWts()
##         p + Upenalty*mean(zth^2)
##     }
## }

## ## Get parameters from deviance function
## ##
## ## @param dfun deviance function created from \code{\link{logisticPcaDevfun}}
## ## @param parms character vectors of parameter names
## ## @param .unlist should the results be \code{\link{unlist}}ed?
## ## @return model parameters at their current value
## ## @export
## getParms <- function(dfun, parms = c("theta", "phi"),
##                      .unlist = TRUE) {
##     ans <- as.list(environment(dfun))[parms]
##     if(.unlist) return(unlist(ans)) else return(ans)
## }

## ## Initialize matrix of column scores
## ##
## ## @param Y response matrix
## ## @param d number of dimensions
## ## @export
## initU <- function(Y, d) {
##     svd2Ym1 <- svd(scale(2*Y - 1))
##     V <- svd2Ym1$v[,1:d, drop = FALSE]
##     qrVt <- qr(t(V))
##     U <- t(qr.R(qrVt))
##     U <- as(U, "dgCMatrix")
##     U@x <- round(U@x, 2)
##     structure(U, sppSVD = V[,1])
## }

## updateInitU <- function(Y1, Y2, U1) {
##     d <- ncol(U1)
##     ## not done obviously
## }

## ## Make logistic PCA model
## ##
## ## @param rho environment of the deviance function
## ## @param opt output of the optimizer
## ## @importFrom Matrix tcrossprod
## ## @export
## mkMod <- function(rho, opt) {
##     s <- rho$s
##     t <- rho$t
##     d <- rho$d
##     theta <- setNames(rho$pp$theta,
##                       c("axis","rows","cols"))
##     beta <- rho$pp$beta(1)
##     b <- rho$pp$b(1)
##     X <- rho$pp$X
##     Zt <- rho$pp$Zt
##     U <- rho$U
##     U@x <- opt$par[-(1:3)]

##     LamU <- Diagonal(t, theta["axis"])
##     LamCol <- Diagonal(t, theta["cols"]^2)
##     typeCors <- tcrossprod(LamU%*%U) + LamCol

##     REindAxes <- 1:(d*s)
##     REindRow <- (d*s+1):((d+1)*s)
##     REindCol <- ((d+1)*s+1):((d+1)*s+t)
    
##     fit <- matrix(rho$pp$linPred(1), s, t)
##     fitInter <- matrix(as.numeric(X%*%beta), s, t)
##     fitAxes <- matrix(as.numeric(b[REindAxes]%*%Zt[REindAxes,]), s, t)
##     fitRow <- matrix(as.numeric(b[REindRow]%*%Zt[REindRow,]), s, t)
##     fitCol <- matrix(as.numeric(b[REindCol]%*%Zt[REindCol,]), s, t)
    
##     svdU <- svd(U, d)
##     colScores <- svdU$u
##     rowScores <- matrix(rho$pp$u(1)[REindAxes], s, d) %*% svdU$v
    
##     return(list(theta = theta, U = U,
##                 typeCors = typeCors,
##                 fit = fit,
##                 fitInter = fitInter,
##                 fitAxes = fitAxes,
##                 fitRow = fitRow,
##                 fitCol = fitCol,
##                 rowScores = rowScores,
##                 colScores = colScores,
##                 beta = beta, b = b,
##                 X = X, Zt = Zt,
##                 REindAxes = REindAxes,
##                 REindRow = REindRow,
##                 REindCol = REindCol))
## }



## ----------------------------------------------------------------------
## glmercModules.R
## ##' Parse mixed model formula for levels-covariance model
## ##'
## ##' @param formula mixed model formula
## ##' @param data data
## ##' @param family family
## ##' @param covList list of covariance matrices for random effects
## ##' grouping factors in \code{formula}
## ##' @param strList list of structure matrices for random effects
## ##' grouping factors in \code{formula}
## ##' @param tileCov in kronecker products, should the covariance
## ##' matrices in \code{covList} be tiled (\code{tileCov = TRUE}) or
## ##' distributed (\code{tileCov = FALSE})?
## ##' @param giveCsparse use compressed-form \code{CsparseMatrix}
## ##' representations (\code{TRUE}), or triplet-form
## ##' \code{TsparseMatrix} representations (\code{FALSE})?
## ##' @param ... not used yet
## ##' @importFrom lme4 findbars nobars
## ##' @export
## glmercFormula <- function(formula, data = NULL, family = binomial,
##                           covList, strList,
##                           tileCov = TRUE,
##                           giveCsparse = TRUE, ...) {

##                                         # get model matrix, grouping
##                                         # factor, and term name
##     reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
##     names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

##                                         # append the covariance and
##                                         # structure matrices over the
##                                         # levels associated with each
##                                         # random effects term
##     for(i in seq_along(reTrmsList)) {
##         reTrmsList[[i]]$strMat <- strList[[reTrmsList[[i]]$grpName]]
##         if(is.null(reTrmsList[[i]]$strMat)) { # if no str mat, use
##                                               # identity
##             nli <- nlevels(reTrmsList[[i]]$grpFac)
##             reTrmsList[[i]]$strMat <- diag(1, nli, nli)
##         }
##         reTrmsList[[i]]$covMat <- covList[[reTrmsList[[i]]$grpName]]
##         if(is.null(reTrmsList[[i]]$covMat)) { # if no cov mat, use
##                                               # identity
##             nstr <- nrow(reTrmsList[[i]]$strMat)
##             reTrmsList[[i]]$covMat <- diag(1, nstr, nstr)
##         }
##     }

##                                         # extract the model matrices,
##                                         # grouping factors, and cov
##                                         # mats
##     modMats <- lapply(reTrmsList, "[[", "modMat")
##     grpFacs <- lapply(reTrmsList, "[[", "grpFac")
##     covMats <- lapply(reTrmsList, "[[", "covMat")
##     strMats <- lapply(reTrmsList, "[[", "strMat")

##                                         # set up constant grouping
##                                         # factors and 1 by 1 cov mats
##                                         # for mkTemplateReTrms (this
##                                         # is where we could add
##                                         # nesting)
##     grpFacConst <- rep(list(as.factor(rep("const", nrow(data)))), length(grpFacs))
##     covMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))
##     strMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))

##                                         # construct the full random
##                                         # effects model matrix,
##                                         # relative covariance factor,
##                                         # etc.
##     if(tileCov) {
##         re <- mkTemplateReTrms(modMats,
##                                grpFacs, grpFacConst,
##                                covMats, covMatConst,
##                                strMats, strMatConst,
##                                giveCsparse)
##     } else {
##         re <- mkTemplateReTrms(modMats,
##                                grpFacConst, grpFacs,
##                                covMatConst, covMats,
##                                strMatConst, strMats,
##                                giveCsparse)
##     }

    
##                                         # organize output value
##     return(c(re,
##              list(X = model.matrix(nobars(formula), data),
##                   y = model.response(model.frame(nobars(formula), data)),
##                   flist = simplifyFacList(grpFacs),
##                   modMats = modMats,
##                   strMats = strMats)))
## }

## ##' Get components of a \code{glmerc} object
## ##'
## ##' @param object \code{\link{glmerc}} object
## ##' @param name object name
## ##' @export
## getMEc <- function(object, name = c("X", "Z", "Zt", "y", "TmodMat", "cnms")) {
##     if (missing(name)) 
##         stop("'name' must not be missing")
##     if (length(name <- as.character(name)) > 1) {
##         names(name) <- name
##         return(lapply(name, getMEc, object = object))
##     }
##     cnms <- lapply(object$parsedForm$modMats, colnames)
##     TmodMat <- object$parsedForm$TmodMat
    
##     for(i in 1:length(TmodMat)) {
##         TmodMat[[i]]@x <- covarByTerms(object)[[i]]
##         dimnames(TmodMat[[i]]) <- rep(cnms[i], 2)
##     }
##     switch(name,
##            X = object$X,
##            Z = t(object$parsedForm$Zt),
##            Zt = object$parsedForm$Zt,
##            y = object$y,
##            TmodMat = TmodMat,
##            cnms = cnms)
## }

## ##' Make covariance template random effects term
## ##'
## ##' @param modMat a model matrix
## ##' @param grpFac1,grpFac2 grouping factors (typically dimension IDs
## ##' of a melted response matrix)
## ##' @param covMat1,covMat2 covariance matrices among the levels of
## ##' \code{grpFac1} and \code{grpFac2}
## ##' @param strMat1,strMat2 structure matrices
## ##' @param giveCsparse use compressed-form \code{CsparseMatrix}
## ##' representations (\code{TRUE}), or triplet-form
## ##' \code{TsparseMatrix} representations (\code{FALSE})?
## ##' @export
## mkTemplateReTrm <- function(modMat,
##                             grpFac1, grpFac2,
##                             covMat1, covMat2,
##                             strMat1, strMat2,
##                             giveCsparse = TRUE) {

##     badDims <- FALSE # FIXME: improve condition
##         # (length(levels(grpFac1)) != nrow(covMat1)) ||
##         # (length(levels(grpFac2)) != nrow(covMat2))
##     if(badDims) stop("covariance matrix incompatible with its grouping factor")

##                                         # use consistent Matrix
##                                         # classes (Triplet-form
##                                         # TSparseMatrix -- i, j, x --
##                                         # not Compressed-form
##                                         # CSparseMatrix)
##     matClass <- if(giveCsparse) "CsparseMatrix" else "TsparseMatrix"

##                                         # expand model matrix from
##                                         # indicator matrices (double
##                                         # `as` required because
##                                         # currently no method in
##                                         # Matrix to go from factor to
##                                         # TsparseMatrix) (FIXME: make
##                                         # sure levels of grpFac's are
##                                         # in the same order as the
##                                         # dimnames of the covMat's)
##     J1 <- as(strMat1 %*% as(grpFac1, "sparseMatrix"), matClass)
##     J2 <- as(strMat2 %*% as(grpFac2, "sparseMatrix"), matClass)
##     Zt <- KhatriRao(KhatriRao(J2, t(modMat)), J1)

##     nc <- ncol(modMat)
##     rowIndices <- sequence(1:nc) # note different order from lme4
##                                  # (FIXME: necessary?)
##     colIndices <- rep(1:nc, 1:nc)
##     TmodMat <- TmodMatLind <- TmodMatOnes <-
##         sparseMatrix(rowIndices, colIndices,
##                      x = 1 * (rowIndices == colIndices),
##                      giveCsparse = giveCsparse)
##     TmodMatLind@x <- as.numeric(seq_along(TmodMat@x))
##     TmodMatOnes@x <- as.numeric(rep(1, length(TmodMat@x)))

##     T1 <- T1ones <- as(chol(covMat1), matClass)
##     T2 <- T2ones <- as(chol(covMat2), matClass)
##     T1ones@x <- as.numeric(rep(1, length(T1ones@x)))
##     T2ones@x <- as.numeric(rep(1, length(T2ones@x)))

##     LambdatBaseline <- T2 %x% TmodMatOnes %x% T1
##     LambdatLind <- T2ones %x% TmodMatLind %x% T1ones
##     LambdatCovar <- T2ones %x% TmodMat %x% T1ones
##     Lambdat <- T2 %x% TmodMat %x% T1

##     return(list(Zt = Zt,
##                 Lambdat = Lambdat, 
##                 LambdatLind = LambdatLind,
##                 LambdatCovar = LambdatCovar,
##                 LambdatBaseline = LambdatBaseline,
##                 TmodMat = TmodMat,
##                 nCovar = length(TmodMat@x)))
## }

## ##' Transpose a list
## ##'
## ##' @param lst \code{\link{list}}
## ##' @export
## listTranspose <- function(lst) {
##     lstExtract <- function(i) lapply(lst, "[[", i)
##     setNames(lapply(names(lst[[1]]), lstExtract), names(lst[[1]]))
## }

## reTrmsBdiag <- function(lst) {
##     ii <- unlist(lapply(lst, slot, "i"))
##     jj <- unlist(lapply(lst, slot, "j"))
##     xx <- unlist(lapply(lst, slot, "x"))
##     sparseMatrix(i = ii, j = jj, x = xx, giveCsparse = FALSE)
## }


## ##' Make several covariance template random effects terms
## ##'
## ##' @param modMat,grpFac1,grpFac2,covMat1,covMat2,strMat1,strMat2,giveCsparse lists of inputs to
## ##' \code{\link{mkTemplateReTrm}}
## ##' @export
## mkTemplateReTrms <- function(modMat,
##                              grpFac1, grpFac2,
##                              covMat1, covMat2,
##                              strMat1, strMat2,
##                              giveCsparse) {

##     reTrmsList <- listTranspose(mapply(mkTemplateReTrm,
##                                        modMat,
##                                        grpFac1, grpFac2,
##                                        covMat1, covMat2,
##                                        strMat1, strMat2,
##                                        SIMPLIFY = FALSE,
##                                        MoreArgs = list(giveCsparse = giveCsparse)))

##     if(!(length(reTrmsList$Zt) == 1L)) {
##         for(i in 2:length(reTrmsList$LambdatLind)) {
##             j <- i - 1
##             reTrmsList$LambdatLind[[i]]@x <-
##                 reTrmsList$LambdatLind[[i]]@x + max(reTrmsList$LambdatLind[[j]]@x)
##         }
##     }
##     within(reTrmsList, {
##         Zt <- do.call(rBind, Zt)
        
##         Lambdat <- .bdiag(Lambdat)
##         LambdatLind <- .bdiag(LambdatLind)
##         LambdatCovar <- .bdiag(LambdatCovar)
##         LambdatBaseline <- .bdiag(LambdatBaseline)
        
##         Lind <- LambdatLind@x
##         covar <- LambdatCovar@x[match(sort(unique(Lind)), Lind)]
##         baseline <- LambdatBaseline@x
##         mapToCovFact <- function(covar) baseline * covar[Lind]
##     })
## }

## ##' Make covariance factor for a covariance template term
## ##'
## ##' @param covMat template covariance matrix
## ##' @param nLevels number of levels over which the grouping factor varies
## ##' @param diagModel templates on the block diagonal (\code{TRUE}) or
## ##' densely tiled over the entire matrix (\code{FALSE})?
## mkTemplateTermLambdat <- function(covMat, nLevels, diagModel = TRUE) {
##     ## cm <- crossprod(matrix(rnorm(25), 5, 5))
##     ## LamtDiag <- mkTemplateTermLambdat(cm, 4, TRUE)
##     ## LamtDense <- mkTemplateTermLambdat(cm, 4, FALSE)
##     ## image(LamtDiag)
##     ## image(LamtDense)
##     Lt <- as(chol(covMat), "sparseMatrix")
##     if(diagModel) {
##         return(.bdiag(rep(list(Lt), nLevels)))
##     } else {
##         n <- nrow(covMat)
##         zeros <- matrix(0, (nLevels - 1) * n, nLevels * n)
##         Lts <- do.call(cBind, rep(list(Lt), nLevels))
##         return(rBind(Lts, zeros))
##     }
## }

## mkTemplateTermZt <- function(explVar, grpFac) {
##     J <- as(as.factor(grp), "sparseMatrix")
##     explVar 
## }

## ##' Make function for computing deviance profiles
## ##' 
## ##' @param mod \code{\link{glmerc}} object
## ##' @param whichPar which parameter to profile
## ##' @return function for computing the profiled deviance over a single
## ##' parameter
## ##' @export
## mkPfun <- function(mod, whichPar) {
##     local({
##         dfun <- mkNfun(mod)
##         optPar <- mod$opt$par
##         whichPar <- whichPar
##         op <- optPar[-whichPar]
##         function(par) {
##             fn <- function(op) {
##                 parNow <- numeric(length(op) + 1)
##                 parNow[whichPar] <- par
##                 parNow[-whichPar] <- op
##                 dfun(parNow)
##             }
##             opt <- minqa:::bobyqa(op, fn, lower = mod$lower[-whichPar],
##                           control = list(iprint = 0L, maxfun = 500))
##             op <<- opt$par # FIXME: assign to local env
##             return(opt$fval)
##         }
##     })
## }

## ##' Make function for computing deviance slices
## ##' 
## ##' @param mod \code{\link{glmerc}} object
## ##' @param whichPar which parameter to slice
## ##' @return function for computing 1D slices through the deviance at
## ##' the optimum
## ##' @export
## mkSfun <- function(mod, whichPar) {
##     local({
##         dfun <- mkNfun(mod)
##         whichPar <- whichPar
##         par <- mod$opt$par
##         function(par1D) {
##             par[whichPar] <- par1D
##             return(dfun(par))
##         }
##     })
## }

## ##' Make function for computing normalized deviance
## ##'
## ##' @param mod model object
## ##' @return function for computing the normalized deviance, which
## ##' equals zero at the optimum
## ##' @export
## mkNfun <- function(mod) {
##     local({
##         dfun <- mod$dfun
##         optFunVal <- dfun(mod$opt$par)
##         function(par) dfun(par) - optFunVal
##     })
## }



## glmerfModules.R
##----------------------------------------------------------------------
## glmerfFormula <- function(formula, data = NULL, family = binomial,
##                           latentName = "latent", latentDims = 0L,
##                           loadingsDim = 1L, loadingPen = 0L,
##                           loadingsPat = NULL, ...) {
    
## }

## ##' Modular functions for generalized bilinear mixed models
## ##'
## ##' @param fr model frame
## ##' @param X fixed effects model matrix
## ##' @param linReTrms linear random effects terms
## ##' @param bilinReTrms bilinear random effects terms
## ##' @param family family
## ##' @param nAGQ nAGQ
## ##' @param verbose verbose
## ##' @param ... ...
## ##' @export
## mkGblmerDevfun <- function(fr, X, linReTrms, bilinReTrms,
##                            family, nAGQ = 1, verbose = 0L, ...) {
##     reTrms <- joinReTrms(bilinReTrms, linReTrms)
##     dfun0  <- mkGlmerDevfun(fr, X, reTrms, family(), nAGQ, verbose, ...)
##     dfun   <- updateGlmerDevfun(dfun0, reTrms)
##     rho    <- environment(dfun)

##     rho$q     <- reTrms$q
##     rho$nth   <- reTrms$nth
##     rho$nCnms <- reTrms$nCnms
##                                         # which Zt@x elements
##                                         # represent loadings
##     rho$Zwhich <- rho$pp$Zt@i %in% (seq_len(rho$q[1]) - 1)
##                                         # mapping from loadings to the
##                                         # Zt@x elements that represent
##                                         # loadings
##     rho$Zind <- with(rho, pp$Zt@x[Zwhich])
##     rho$Ztx <- rho$pp$Zt@x
    
##     return(dfun)
## }

##----------------------------------------------------------------------
## glmerc.R

## ##' Generalized linear mixed model with custom covariance over grouping factor levels
## ##'
## ##' When fitting generalized linear mixed model with
## ##' \code{\link{glmer}}, the covariance structure over the grouping
## ##' factor levels is assumed to be an identity matrix.  That is,
## ##' random effects are sampled indepedently over the grouping factor
## ##' levels.  With \code{glmerc}, one may specify the covariance over
## ##' the levels, up to a fitted parameter.
## ##'
## ##' @param formula formula
## ##' @param data data
## ##' @param family family
## ##' @param covList list of covariance matrices tagged by the names of
## ##' grouping factors
## ##' @param strList list of structure matrices tagged by the names of
## ##' grouping factors
## ##' @param tileCov in kronecker products, should the covariance
## ##' matrices in \code{covList} be tiled (\code{tileCov = TRUE}) or
## ##' distributed (\code{tileCov = FALSE})?
## ##' @param giveCsparse use compressed-form \code{CsparseMatrix}
## ##' representations (\code{TRUE}), or triplet-form
## ##' \code{TsparseMatrix} representations (\code{FALSE})?
## ##' @param optControl optControl
## ##' @param ... ...
## ##' @export
## glmerc <- function(formula, data = NULL, family = binomial,
##                    covList = list(), strList = list(),
##                    tileCov = TRUE,
##                    giveCsparse = TRUE,
##                    optControl = list(iprint = 0L), ...) {

##                                         # parse formula
##     data <- as.data.frame(data)
##     parsedForm <- glmercFormula(formula, data,
##                                 covList = covList,
##                                 strList = strList,
##                                 tileCov = tileCov,
##                                 giveCsparse = giveCsparse)

##                                         # organize initial values
##     covar <- parsedForm$covar
##     fixef <- rep(0, ncol(parsedForm$X))
##     initPars <- c(covar = covar,
##                   fixef = fixef)
##     parInds <- list(covar = seq_along(covar),
##                     fixef = seq_along(fixef) + length(covar),
##                     loads = NULL)

##                                         # construct deviance function
##     dfun <- mkGeneralGlmerDevfun(parsedForm$y, parsedForm$X,
##                                  parsedForm$Zt, parsedForm$Lambdat,
##                                  rep(1, nrow(data)), rep(0, nrow(data)),
##                                  initPars, parInds,
##                                  parsedForm$mapToCovFact, function(loads) NULL,
##                                  ...)

##                                         # optimize deviance function
##     dfun(initPars)
##     lower <- ifelse(initPars, 0, -Inf)
##     opt <- minqa:::bobyqa(initPars, dfun, lower = lower,
##                   control = optControl)
##     if(FALSE) {for(i in 1:5) {
##         opt <- minqa:::bobyqa(opt$par, dfun, lower = lower,
##                       control = optControl)
##     }}
##     names(opt$par) <- names(initPars)

##                                         # organize return value
##     ans <- list(opt = opt, parsedForm = parsedForm, dfun = dfun,
##                 parInds = parInds, X = parsedForm$X, y = parsedForm$y,
##                 lower = lower)
##     class(ans) <- "glmerc"
##     return(ans)
## }

## ##' @param ... not used
## ##' @importFrom lme4 fixef
## ##' @rdname pars
## ##' @export
## fixef.glmerc <- function(object, ...) {
##     setNames(.fixef(object$opt$par, object$parInds),
##              colnames(object$X))
## }

## ##' @rdname pars
## ##' @export
## covar.glmerc <- function(object, ...) .covar(object$opt$par, object$parInds)


## ##' @rdname pars
## ##' @export
## loads.glmerc <- function(object, ...) stop("covariance over levels models do not have loadings")


## ##' @rdname pars
## ##' @export
## pars.glmerc <- function(object, ...) c(covar(object), fixef(object))

## ##' @rdname pars
## ##' @export
## covarByTerms <- function(object, ...) {
##     nCovar <- object$parsedForm$nCovar
##     if(length(nCovar) == 1L) return(list(covar(object)))
##     inds <- covarInds(object$parsedForm$nCovar)
##     lapply(inds, function(i) covar(object)[i])
## }

## covarInds <- function(nCovar) {
##     ans <- lapply(nCovar, seq, from = 1, by = 1)
##     for(i in 2:length(nCovar)) ans[[i]] <- ans[[i]] + max(ans[[i-1]])
##     return(setNames(ans, names(nCovar)))
## }

## ##' @importFrom lme4 VarCorr
## ##' @export
## VarCorr.glmerc <- function(x, ...) {
##     tmm <- getMEc(x, "TmodMat")
##     cnms <- getMEc(x, "cnms")
##     lapply(lapply(tmm, as.matrix), crossprod)
## }

## ##' @export
## vcov.glmerc <- function(object, justFixef = TRUE, ...) {
##     ## FIXME: DRY with vcov.strucGlmer
##     optPar <- object$opt$par
##     ans <- solve(0.5 * lme4:::deriv12(object$dfun,
##                                       optPar)$Hessian)
##     dimnames(ans) <- rep(list(names(optPar)), 2)
##     if(justFixef) {
##         dims <- object$parInds$fixef
##         ans <- ans[dims, dims]
##     }
##     return(ans)
## }

## .safeExtractDiag <- function(x) {
##     if(length(x) == 1) return(x)
##     diag(x)
## }

## ##' @param x \code{glmerc} object
## ##' @rdname glmerc
## ##' @export
## print.glmerc <- function(x, ...) {
##     cat("\nGeneralized linear mixed model\nwith covariance amongst grouping factor levels\n")
##     cat("----------------------------------------------\n\n")

##     cat("Fixed effects\n")
##     cat("-------------\n\n")
##     print(cbind(Estimate = fixef(x),
##                 `Std. Error` = sqrt(.safeExtractDiag(vcov(x)))))

##     cat("\n\nRandom effects (co)variance\n")
##     cat("---------------------------\n\n")
##     print(VarCorr.glmerc(x))
## }


## ----------------------------------------------------------------------
## glmerf.R
## ##' Generalized linear mixed model with factor loadings
## ##'
## ##' @param formula a mixed model formula for the linear part of the
## ##' model
## ##' @param data a \code{\link{data.list}} object
## ##' @param family a \code{\link{family}} object
## ##' @param latentName name for latent factor
## ##' @param latentDims number of latent dimensions in the bilinear
## ##' component of the model
## ##' @param loadingsDim what dimension of \code{data} is associated
## ##' with loadings?
## ##' @param loadingPen penalty for the size of the loadings
## ##' @param verbose passed to optimizer
## ##' @param control arguments for the optimizer
## ##' @param ... arguments to be passed to \code{\link{glFormula}}
## ##' @importClassesFrom lme4 merMod glmerMod
## ##' @import multitable
## ##' @importFrom lme4 glmer glFormula mkGlmerDevfun updateGlmerDevfun
## ##' @export
## glmerf <- function(formula, data, family,
##                    latentName = "latent", latentDims = 0L,
##                    loadingsDim = 1L, loadingPen = 0L,
##                    verbose = 0L, control = list(), ...) {
##     if(!any(inherits(data, "data.list"))) stop("data must be a data list")
##     if(length(dim(data)) != 2L) stop("data list must be two dimensional")
##     dIds <- names(dd <- dim(data))
##     df <- as.data.frame(dims_to_vars(data))
##     initGlmer <- glmer(formula, df, family, ...)
##     if(latentDims == 0L) {
##         warning("no latent variables, returning glmer results")
##         return(initGlmer)
##     }
##     if(loadingsDim == 2L) {
##         datY <- data$Y
##     } else if(loadingsDim == 1L) {
##         datY <- t(data$Y)
##     } else {
##         stop("loadingsDim neight 1 nor 2")
##     }
##     U <- try(matrix(0, nrow = dd[loadingsDim], ncol = latentDims))
##     if(inherits(U, "try-error")) stop("loadingsDim does not index a dimension of data")
##     nFreeLoadings <- (dd[loadingsDim] * latentDims) - choose(latentDims, 2)
##     U[lower.tri(U, TRUE)] <- 1:nFreeLoadings
##     latentVarNames <- paste(latentName, 1:latentDims, sep = "")
##     U <- setNames(as.data.frame(U), latentVarNames)
##     latentData <- data.list(U, drop = FALSE, dimids = dIds[loadingsDim])
##     data <- data + latentData
##     df <- as.data.frame(dims_to_vars(data))

##     # form1 <- formula
##     # form2 <- as.formula(paste(". ~ 0 + (0 + latent |", names(dd[2], ")"))
##     linFormula <- formula
##     latentForm <- paste(latentVarNames, collapse = " + ")
##     bilinFormula <- as.formula(paste(". ~ 0 + (0 + ", latentForm, " ||",
##                                      dIds[-loadingsDim], ")"))
##                                         # FIXME: replace -tv with more
##                                         # general removal strategy
##                                         # (e.g. with characters etc)
##     bilinFormula[[2]] <- linFormula[[2]] # use same response variables

##       parsedLinFormula <- parsedForm <- glFormula(  linFormula, df, family, ...)
##     parsedBilinFormula <-               glFormula(bilinFormula, df, family, ...)
##     reTrms <- joinReTrms(parsedBilinFormula$reTrms, parsedLinFormula$reTrms)
##     parsedForm$reTrms <- reTrms

##     theta <- reTrms$theta
##     lower <- reTrms$lower
##     dfun <- do.call(mkGlmerDevfun, parsedForm)
##     dfun <- updateGlmerDevfun(dfun, reTrms)
##     rho <- environment(dfun)

##                                         # which Zt@x elements
##                                         # represent loadings
##     rho$Zwhich <- rho$pp$Zt@i %in% (seq_len(latentDims * dd[-loadingsDim]) - 1)
##                                         # mapping from loadings to the
##                                         # Zt@x elements that represent
##                                         # loadings
##     rho$Zind <- with(rho, pp$Zt@x[Zwhich])
##     rho$Ztx <- rho$pp$Zt@x
##     rho$loadInd <- 1:nFreeLoadings

##     dfunPrefix <- function(pars) {
##         loadings <- pars[loadInd]
##         Ztx[Zwhich] <- loadings[Zind]
##         trash <- pp$setZt(Ztx)
##         pars <- c(rep(1, latentDims), pars[-loadInd])
##                                         # the '1' is is for scalar
##                                         # bilinear random effects
##                                         # (FIXME: be more general)
##     }

##     dfunSuffix <- function(pars) {
##         p + sum(loadingPen * abs(loadings))
##     }

##     body(dfun) <- cBody(body(dfunPrefix), body(dfun), body(dfunSuffix))
##     formals(dfun) <- setNames(formals(dfun), "pars")

##     initLoadings <- svd(scale(datY))$v[, 1:latentDims, drop = FALSE]
##     initLoadings <- initLoadings[lower.tri(initLoadings, TRUE)]
    
##     #opt <- optim(c(initLoadings, theta[-1]), dfun, method = "L-BFGS-B",
##     #             lower = c(rep(-Inf, dd[1]), lower),
##     #             control = list(trace = 3))

##                                         # here is the '1' again (with
##                                         # a minus in front) for scale
##                                         # bilinear random effects

##     # initial parameters
##     # order: (1) loadings, (2) covariance pars, (3) fixed effect pars
##     parPointer <- setNames(cumsum(c(0, length(initLoadings),
##                                     length(theta) - latentDims)) + 1,
##                            c("loadings", "theta", "fixef"))
##     initPars <- c(initLoadings, theta[-(1:latentDims)], rho$pp$beta(1))
##     optLower <- c(rep(-Inf, length(initLoadings)), rho$lower[-(1:latentDims)])
##                                         # optimize
##     opt <- lme4:::optwrap("bobyqa", dfun, initPars, 
##                           lower = optLower, verbose = verbose,
##                           control = control)

##     optPar <- opt$par
##     optNoLoadings <- c(rep(1, latentDims), opt$par[-rho$loadInd])
##     optLoadings <- matrix(0, dd[loadingsDim], latentDims)
##     optLoadings[lower.tri(optLoadings, TRUE)] <- optPar[rho$loadInd]
##     dim(optLoadings) <- c(dd[loadingsDim], latentDims)
##     colnames(optLoadings) <- latentVarNames
##     rownames(optLoadings) <- dimnames(data)[[loadingsDim]]
##     opt$par <- optNoLoadings    

##     ## rho$control <- attr(opt,"control")
##     # rho$nAGQ <- 0

##     # body(dfun) <- body(dfun)[c(1, 5:9)]
##     # formals(dfun) <- setNames(formals(dfun), "theta")

##     ## opt <- optimizeGlmer(dfun) # optimize without updating loadings
##     ## opt$par <- optTheta
##     # dfun(optTheta)

##     mer <- mkMerMod(environment(dfun), opt, parsedForm$reTrms, parsedForm$fr)

##                                         # (FIXME: write specific
##                                         # mkGlmerLatentMod)
##     merList <- list()
##     for(sl in names(getSlots("glmerMod"))) {
##         merList[[sl]] <- slot(mer, sl)
##     }
##     do.call(new, c(list(Class = "gblmerMod"),
##                    merList,
##                    list(  loadings = optLoadings,
##                             optPar = optPar,
##                         parPointer = parPointer)), quote = TRUE)
## }

## ##' Class "gblmerMod"
## ##'
## ##' Class \code{"gblmerMod"} is san S4 class that extends
## ##' \code{"glmerMod"}
## ##' @name gblmerMod-class
## ##' @aliases gblmerMod
## ##' @aliases gblmerMod-class
## ##' \section{Slots}{
## ##'   \describe{
## ##'     \item{\code{loadings}:}{factor loadings.}
## ##'   }
## ##' }
## ##' @keywords classes
## ##' @export
## setClass("gblmerMod",
##          representation(  loadings = "matrix",
##                             optPar = "numeric",
##                         parPointer = "numeric"),
##          contains="glmerMod")

## ##' Variance-covariance matrix of the (co)variance parameters and
## ##' loadings (and sometimes of fixed effects)
## ##'
## ##' @param object a \code{\link{gblmerMod}} object
## ##' @param correlation should the correlation matrix be computed too?
## ##' @param ... not used
## ##' @export
## vcov.gblmerMod <- function(object, correlation = TRUE, ...) {
##                                         # calc.vcov.hess is a function
##                                         # from vcov.merMod (FIXME:
##                                         # breakout of vcov.merMod and
##                                         # expose?)
##     calc.vcov.hess <- function(h) {
## 	## ~= forceSymmetric(solve(h/2)[i,i]) : solve(h/2) = 2*solve(h)
##         h <- tryCatch(solve(h),
##                       error=function(e) matrix(NA,nrow=nrow(h),ncol=ncol(h)))
##         ## i <- -seq_len(ntheta)
## 	## h <- h[i,i]
## 	forceSymmetric(h + t(h))
##     }
##     fixefIndices <- object@parPointer[3]:length(object@optPar)
##     h <- object@optinfo$derivs$Hessian[fixefIndices, fixefIndices, drop = FALSE]
##     V.hess <- calc.vcov.hess(h)
    
##     bad.V.hess <- any(is.na(V.hess))
##     if (!bad.V.hess) {
##         e.hess <- eigen(V.hess,symmetric = TRUE,only.values = TRUE)$values
##         if (min(e.hess) <= 0) bad.V.hess <- TRUE
##     }
##     if (!bad.V.hess) {
##         V <- V.hess
##     } else {
##         stop("variance-covariance matrix computed ",
##              "from finite-difference Hessian is\n",
##              "not positive definite or contains NA values")
##     }
##     rr <- tryCatch(as(V, "dpoMatrix"), error = function(e)e)
##     if (inherits(rr, "error")) {
## 	warning(gettextf("Computed variance-covariance matrix problem: %s;\nreturning NA matrix",
##                          rr$message), domain = NA)
##         rr <- matrix(NA,nrow(V),ncol(V))
##     }
##     if(correlation)
## 	rr@factors$correlation <-
## 	    as(rr, "corMatrix") else rr # (is NA anyway)
##     rr
## }


## ##' Refactor \code{loadings} into a generic function
## ##'
## ##' @param x object with loadings to extract
## ##' @param ... additional arguments
## ##' @export
## loadings <- function(x, ...) {
##     UseMethod("loadings")
## }

## ##' @export
## loadings.default <- function(x, ...) x$loadings

## ##' @export
## loadings.gblmerMod <- function(x, ...) x@loadings

## ##' Concatenate the bodies of functions
## ##'
## ##' @param ... function bodies to combine
## ##' @export
## cBody <- function(...) {
##     l... <- list(...)
##     l...[-1] <- lapply(l...[-1], "[", -1)
##     l...asList <- lapply(l..., as.list)
##     l...cList <- do.call(c, l...asList, quote = TRUE)
##     as.call(l...cList)
## }

## ##' Join two lists describing two sets of random effects terms
## ##'
## ##' @param reTrms1,reTrms2 results of \code{\link{glFormula}(.)$reTrms}
## ##' @export
## joinReTrms <- function(reTrms1, reTrms2) {
##     names(reTrms1) <- paste(names(reTrms1), 1, sep = "")
##     names(reTrms2) <- paste(names(reTrms2), 2, sep = "")
##     reTrms <- c(reTrms1, reTrms2)
##     with(reTrms, {
##         flist <- joinFlist(flist1, flist2)
##         attr(flist, "assign") <- c(match(names(cnms1), names(flist)),
##                                    match(names(cnms2), names(flist)))
##         q <- c(nrow(Zt1), nrow(Zt2))
##         nth <- c(length(theta1), length(theta2))
##         nCnms <- c(length(cnms1), length(cnms2))
##         list(Zt = rBind(Zt1, Zt2),
##              theta = c(theta1, theta2),
##              Lind = c(Lind1, Lind2 + max(Lind1)),
##              Gp = c(Gp1, Gp2[-1] + max(Gp1)),
##              lower = c(lower1, lower2),
##              Lambdat = .bdiag(list(Lambdat1, Lambdat2)),
##              flist = flist,
##              cnms = c(cnms1, cnms2),
##              Ztlist = c(Ztlist1, Ztlist2),
##              q = q, nth = nth, nCnms = nCnms)
##     })
## }

## ## ------------------------------------------------------------
## ## adpated from mkReTrms
## ## ------------------------------------------------------------
## joinFlist <- function(flist1, flist2) {
##     flist <- c(flist1, flist2)
##     fnms <- names(flist)
##     if (length(fnms) > length(ufn <- unique(fnms))) {
##         flist <- flist[match(ufn, fnms)]
##         asgn <- match(fnms, ufn)
##     } else asgn <- seq_along(flist)
##     names(flist) <- ufn
##     flist <- do.call(data.frame, c(flist, check.names = FALSE))
##     return(flist)
## }

## ##' Get random effects terms from a fitted \code{merMod} object
## ##'
## ##' @param object \code{\link{merMod}} object
## ##' @param ... not used
## ##' @return see \code{\link{mkReTrms}}
## ##' @importFrom lme4 getME
## ##' @export
## getReTrms <- function(object, ...) {
##     rts <- c("Zt", "Lambdat", "Lind", "theta", "lower", "flist", "cnms")
##     getME(object, rts)
## }

