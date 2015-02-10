library(Matrix)
library(lme4ord)
library(lme4)
library(lme4pureR)
library(multitable)
library(pryr)
library(reo)

## n <- 50
## m <- 30
## x <- rnorm(n)
## z <- rnorm(m)
## B <- rmat(2, 2)
## u <- rnorm(n)
## v <- rnorm(m)

## eta <-
##     (cbind(1, x) %*% B %*% rbind(1, z)) +
##     (cbind(1, x) %*%        rmat(2, m)) +
##     ( rmat(n, 2) %*%       rbind(1, z)) +
##     (   cbind(u) %*%          rbind(v))

## Y <- rmat(n, m,
##           f(rbinom(n, size = 1, prob = plogis(eta))),
##           eta = eta)

data(fish)
data(limn)
Y <- as.matrix(fish)
## Y <- Y[, colSums(Y) > 1]
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))
dl <- data.list(Y = Y, x = x,
                dimids = c("sites", "species"))
summary(dl)

mod <- gblmer(Y ~ 1 + (1 | species), dl, binomial,
              loadingsDim = 2,
              latentDim = 2,
              loadingPen = 1L,
              verbose = 2L)

## dl <- dl + variable(1:dim(dl)[1], "species", "latent")
## df <- as.data.frame(dims_to_vars(dl))
##   linFormula <- Y ~ 1 + (1          | species)
## bilinFormula <- Y ~ 0 + (0 + latent | sites  )
##   parsedLinFormula <- glFormula(  linFormula, df, binomial)
## parsedBilinFormula <- glFormula(bilinFormula, df, binomial)
## reTrms <- joinReTrms(parsedBilinFormula$reTrms, parsedLinFormula$reTrms)
## dfun <- mkGblmerDevfun(parsedLinFormula$fr, parsedLinFormula$X,
##                          parsedLinFormula$reTrms,
##                        parsedBilinFormula$reTrms, binomial)
## rho <- environment(dfun)
## rho$pp$Zt@i %in% (seq_len(rho$q[1]) - 1)

## mod <- gblmer(Y ~ 1 + (1 | species), . ~ 0 + (0 + latent | sites),
##               dl, binomial, 1, 1, 2)

mod
summary(mod)

ranef(mod)$species
ranef(mod)$sites
(loadMod <- loadings(mod))
latentCov <- Matrix(loadMod %*% t(loadMod)) + diag(VarCorr(mod)$species[1], m, m)
image(cov2cor(latentCov))
image(vcov(mod)@factors$correlation)


fitY <- matrix(getME(mod, "mu"), n, m, byrow = TRUE)
boxplot(fitY ~ Y, las = 1, horizontal = TRUE)

## j <- 3
## plot(fitY[,j], Y[,j])

# plot(u, getME(mod, "b")[1:n])
# plot(v, mod@loadings)

j <- 14
b <- getME(mod, "b")[1:n]
bb <- seq(min(b), max(b), length = 100)
plot(b, Y[,j])
lines(bb, plogis(loadings(mod)[j]*bb))


