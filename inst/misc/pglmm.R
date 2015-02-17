library(pez)
library(lme4ord)
library(minqa)
library(plotrix)

n <- 10
m <- 30
dl <- dims_to_vars(data.list(y = 1 * (rmat(n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)
covList <- list(species = stanCov(cov(matrix(rnorm(m * (m + 1)), m + 1, m))))

parsedForm <- levelsCovFormula(y ~ x + (x | species), df, covList = covList)

parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(c(0.5, -0.2, 0.5)))

X <- model.matrix(y ~ x, df)
Z <- t(parsedForm$Lambdat %*% parsedForm$Zt)
beta <- rnorm(ncol(X))
u <- rnorm(ncol(Z))
p <- plogis(as.numeric(X %*% beta + Z %*% u))
dl$y <- rbinom(nrow(df), 1, p)

hist(p)

parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(c(1, 0, 1)))

color2D.matplot(dl$y, xlab = "species", ylab = "sites", main = "abundance")

image(parsedForm$Lambdat)

df <- as.data.frame(dl)

parInds <- list(covar = 1:3, fixef = 4:5, loads = NULL)
initPars <- c(covar = c(1, 0, 1), fixef = c(0, 0))
dfun <- mkGeneralGlmerDevfun(parsedForm$y, parsedForm$X,
                             parsedForm$Zt, parsedForm$Lambdat,
                             rep(1, nrow(df)), rep(0, nrow(df)),
                             initPars, parInds,
                             parsedForm$mapToCovFact, function(loads) NULL)
dfun(initPars)

opt <- bobyqa(initPars, dfun, lower = c(0, -Inf, 0, -Inf, -Inf),
              control = list(iprint = 4L))
names(opt$par) <- names(initPars)
opt$par








modMat  <- list(model.matrix( ~ 0 + x, df),
                model.matrix( ~ 1, df))
grpFac1 <- list(df$sites, df$species)
grpFac2 <- list(df$const, df$const)
covMat1 <- list(crossprod(rmat(n + 1, n)),
                crossprod(rmat(m + 1, m)))
covMat2 <- list(matrix(1, 1, 1),
                matrix(1, 1, 1))

ret <- mkTemplateReTrms(list(model.matrix( ~ x, df)),
                        list(df$species),
                        list(df$const),
                        list(crossprod(rmat(m + 1, m))),
                        list(matrix(1, 1, 1)))

for(i in 1:3) dev.new()

#ret <- mkTemplateReTrms(modMat, grpFac1, grpFac2, covMat1, covMat2)
#ret <- within(ret, Lambdat@x[] <- mapToCovFact(c(1, 0.5)))

ret <- within(ret, Lambdat@x[] <- mapToCovFact(c(0.1, -0.05, 0.1)))

X <- model.matrix(Y ~ x, df)
Z <- t(ret$Lambdat %*% ret$Zt)
beta <- rnorm(ncol(X))
u <- rnorm(ncol(Z))
p <- plogis(as.numeric(X %*% beta + Z %*% u))
dl$Y <- rbinom(nrow(df), 1, p)

hist(p)

ret <- within(ret, Lambdat@x[] <- mapToCovFact(c(1, 0, 1)))

image(dl$Y)
color2D.matplot(dl$Y, xlab = "species", ylab = "sites", main = "abundance")

image(ret$Lambdat)

df <- as.data.frame(dl)

## dev.set(2)
## image(ret$Lambdat)
## image(environment(dfun)$pp$Lambdat)
## #image(ret$Lambdat[1:12, 1:12])
## #image(ret$Lambdat[13:18, 13:18])
## dev.set(3)
## image(ret$Zt)
## dev.set(4)
## image(crossprod(ret$Lambdat %*% ret$Zt))

parInds <- list(covar = 1:3, fixef = 4:5, loads = NULL)
initPars <- c(covar = c(1, 0, 1), fixef = c(0, 0))
dfun <- mkGeneralGlmerDevfun(df$Y, X,
                             ret$Zt, ret$Lambdat,
                             rep(1, nrow(df)), rep(0, nrow(df)),
                             initPars, parInds,
                             ret$mapToCovFact, function(loads) NULL)
dfun(initPars)

opt <- bobyqa(initPars, dfun, lower = c(0, -Inf, 0, -Inf, -Inf),
              control = list(iprint = 4L))
names(opt$par) <- names(initPars)

dfun(opt$par)
opt$par
rho <- environment(dfun)
image(rho$pp$Lambdat)
image(ret$Lambdat)

ret$Lambdat@x[] <- ret$mapToCovFact(opt$par)

covEta <- t(ret$Zt) %*% crossprod(ret$Lambdat) %*% ret$Zt
image(cov2cor(covEta))
image(covEta)
image(ret$Lambdat)
image(crossprod(ret$Lambdat))
