library(lme4ord)

############################################################
## silly simulation example
############################################################

n <- 60 # samples 
p <- 2  # fixed effects
q <- 30 # random effects

                                        # covariance factor
covTemplate <- as(chol(cov(matrix(rnorm((q+1)*q), q + 1, q))),
                  "sparseMatrix")
Lambdat <- covTemplate
                                        # random effects model matrix
Zt <- sparseMatrix(i = rep(1:q, 6),
                   j = as.vector(outer(rep(1:10, 3), seq(0, 50, 10), "+")),
                   x = rep(rnorm(10), 18))
                                        # fixed effects model matrix
X <- matrix(rnorm(n * p), n, p)
                                        # response vector
oeta <- as.numeric(X %*% c(1, 1) + t(Zt) %*% t(Lambdat) %*% rnorm(q))
y <- rbinom(n, 1, plogis(eta))
weights <- rep(1, n)
offset <- rep(0, n)


initPars <- c(covar = 1,
              fixef = c(0, 0),
              loads = rnorm(10))
parInds <- list(covar = 1,
                fixef = 2:3,
                loads = 4:13)

                                        # mappings for updating
                                        # parameters
mapToCovFact <- function(covar) covar * covTemplate@x
mapToModMat <- function(loads) rep(loads, 18)
                                        # e.g. set to initial values
Lambdat@x <- mapToCovFact(initPars[parInds$covar])
Zt@x <- mapToModMat(initPars[parInds$loads])

image(Zt)
image(Lambdat)

devfun <- mkGeneralGlmerDevfun(y, X, Zt, Lambdat,
                               weights, offset,
                               initPars, parInds,
                               mapToCovFact, mapToModMat)                               
                               

devfun(initPars)

opt <- minqa:::bobyqa(initPars, devfun)
setNames(opt$par, names(initPars))
