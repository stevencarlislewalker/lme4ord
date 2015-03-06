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
