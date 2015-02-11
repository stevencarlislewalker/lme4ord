##' Make random effects structures for the single correlation template model
##'
##' @param corr a correlation matrix template
##' @param grp a grouping factor vector
##' @param n sample size
##' @export

mkLambdatTemplateModel <- function(covMat, nLevels, blockDiag = TRUE, upper = NULL, lower = NULL) {
    if(missing(covMat) & (!is.null(upper))) {
        covFact <- t(upper)
    } else if(missing(covMat) & (!is.null(lower))) {
        covFact <- lower
    } else {
        covFact <- t(chol(covMat))
    }
    structMat <- if(blockDiag) diag(1, nLevels, nLevels) else matrix(1, nLevels, nLevels)
    as(kronecker(structMat, covMat), Class = "sparseMatrix")
}

##image(mkLambdatTemplateModel(re.2$covar, 10))
cc <- mkLambdatTemplateModel(re.2$covar, 2, FALSE)
##image(cc)
##Cholesky(cc, TRUE, LDL = TRUE, super = FALSE)
eigcc <- eigen(cc)
ff <- eigcc$vectors %*% diag(sqrt(round(eigcc$values, 4)))
ffqr <- qr(t(ff))
##image(as(round(t(qr.R(ffqr)) %*% qr.R(ffqr), 10), "sparseMatrix"))
##image(cc)
image(as(round(qr.R(ffqr), 10), "sparseMatrix"))
round(qr.R(ffqr)[1:15, 16:30], 2)
round(qr.R(ffqr)[1:15, 1:15], 2)
image(t(qr.R(ffqr)))

ffR <- as(round(t(qr.R(ffqr)), 10), "sparseMatrix")
image(t(ffR) %*% ffR)

      
eigcc0 <- eigen(cc[1:15, 1:15])
ff0 <- as(eigcc0$vectors %*% diag(eigcc0$values), "sparseMatrix")
image(ff0)
image(ff[1:15, 1:15])
image(as(ff %*% t(ff), "sparseMatrix"))
image(as(cc, "sparseMatrix"))
qreigcc <- qr(t(eigcc$vectors))
fac <- t(qr.R(qreigcc)) %*% t(qr.Q(qreigcc)) %*% diag(sqrt(round(eigcc$values, 5)))
image(fac)

 %*% qr.Q(qreigcc) %*% qr.R(qreigcc)
eigcc$vectors

svdMLT$u[,1:15] %*% diag(svdMLT$d[1:15])


image(chol(1 * mkLambdatTemplateModel(re.2$covar, 10, TRUE) +
    runif(1) * mkLambdatTemplateModel(re.2$covar, 10, FALSE)))

mkZtTemplateModel <- function(modelMatrix, groupFac, blockDoag = TRUE) {
    
}



XX <- t(XX) %*% XX
kronecker(matrix(1, 


image(as(kronecker(matrix(1, 10, 15), t(chol(re.2$covar))), Class = "sparseMatrix"))


mkRanefStructuresCorr <- function(corr, grp, n){
                                        # create indicator matrix and order it
                                        # to be consistent with the order of corr
    Jt <- as(as.factor(grp), Class="sparseMatrix")
    Jt <- Jt[dimnames(corr)[[1]],]

                                        # create Zt
    Zt <- KhatriRao(Jt, t(rep(1,n)))

                                        # create Lambdat
    Lambdat <- as(t(chol(corr)), Class="sparseMatrix")

                                        # create mapping from theta to
                                        # the non-zero components of
                                        # Lambdat
    thfun <- local({
        template <- Lambdat
        function(theta) theta * template@x})

    list(     Zt = Zt,
         Lambdat = Lambdat,
           thfun = thfun,
           theta = 1,
           lower = 0,    # lower and
           upper = Inf)  # upper bounds on theta parameters
}
