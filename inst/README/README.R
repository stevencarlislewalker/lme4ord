
## ----, setup, echo = FALSE-----------------------------------------------
library(knitr)
opts_chunk$set(fig.path = "inst/README/figure/")


## ----, packages, message = FALSE-----------------------------------------
library(Matrix)
library(lme4ord)
library(plotrix)
library(minqa)
library(ape)
library(lme4)
library(lme4pureR)
library(multitable)
library(pryr)
library(reo)


## ----, fig.width = 4, fig.height = 6-------------------------------------
td <- simTestPhyloDat(10, n = 20, m = 10, power = 0.4)
color2D.matplot(td$dl$y, xlab = "species", ylab = "sites",
                main = "Occurrence")
plot(td$ph)
edgelabels()


## ------------------------------------------------------------------------
(indMat <- edgeTipIndicator(td$ph))


## ------------------------------------------------------------------------
dummy <- as.data.frame(t(indMat))
td$dl <- td$dl + variableGroup(dummy, "species")
edgeNms <- names(dummy)
df <- as.data.frame(td$dl)


## ------------------------------------------------------------------------
Z <- model.matrix(as.formula(paste("~ 0 + ", paste(edgeNms, collapse = " + "))), df)
X <- model.matrix(~ 1, df)
y <- model.response(model.frame(y ~ 1, df))
n <- nrow(df)
p <- ncol(X)
q <- ncol(Z)
mapToCovFact <- local({
    q <- q
    function(covar) rep(covar, q)
})


## ------------------------------------------------------------------------
dfun <- mkGeneralGlmerDevfun(y = y, X = X,
                             Zt = as(t(Z), "sparseMatrix"),
                             Lambdat = sparseMatrix(i = 1:q, j = 1:q, x = 1),
                             weights = rep(1, n), offset = rep(0, n),
                             initPars = c(1, 0),
                             parInds = list(covar = 1, fixef = 2),
                             mapToCovFact = mapToCovFact,
                             mapToModMat = NULL)


## ------------------------------------------------------------------------
opt <- optim(c(1, 0), dfun, lower = c(0, -Inf), method = "L-BFGS-B")
dfun(opt$par)
opt$par


## ----, fig.height = 4, fig.width = 4-------------------------------------
rho <- environment(dfun)
with(rho$pp, image(crossprod(Lambdat %*% Zt)))


## ----, fig.height = 8, fig.width = 4-------------------------------------
plot(td$ph)
edgelabels(round(rho$pp$b(1), 2), cex = 1)


## ----, echo = 15---------------------------------------------------------
td <- simTestPhyloDat(10, n = 100, m = 500, power = 0.1)
indMat <- edgeTipIndicator(td$ph)
dummy <- as.data.frame(t(indMat))
td$dl <- td$dl + variableGroup(dummy, "species")
edgeNms <- names(dummy)
df <- as.data.frame(td$dl)
Z <- model.matrix(as.formula(paste("~ 0 + ", paste(edgeNms, collapse = " + "))), df)
X <- model.matrix(~ 1, df)
y <- model.response(model.frame(y ~ 1, df))
n <- nrow(df)
p <- ncol(X)
q <- ncol(Z)
mapToCovFact <- local({
    q <- q
    function(covar) rep(covar, q)
})
dfun <- mkGeneralGlmerDevfun(y = y, X = X,
                             Zt = as(t(Z), "sparseMatrix"),
                             Lambdat = sparseMatrix(i = 1:q, j = 1:q, x = 1),
                             weights = rep(1, n), offset = rep(0, n),
                             initPars = c(1, 0),
                             parInds = list(covar = 1, fixef = 2),
                             mapToCovFact = mapToCovFact,
                             mapToModMat = NULL)
system.time(opt <- optim(c(1, 0), dfun, lower = c(0, -Inf), method = "L-BFGS-B"))


## ----, fig.width = 4, fig.height = 4-------------------------------------
image(as(tcrossprod(indMat), "sparseMatrix"))


## ----, initial sims------------------------------------------------------
set.seed(10)
n <- 10
m <- 30
dl <- dims_to_vars(data.list(y = 1 * (matrix(rnorm(n * m), n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)


## ----, phylogeny sims----------------------------------------------------
phy <- rtree(n = m)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)


## ----, covariance sims---------------------------------------------------
Vphy <- stanCov(vcv(phy))
dimnames(Vphy) <- rep(list(1:m), 2)


## ----, plot phylogeny, fig.width = 6.5, fig.height = 6.5-----------------
plot(phy)
image(as(Vphy, "sparseMatrix"))


## ----, make the covList--------------------------------------------------
covList <- list(species = Vphy)


## ----, formula parsing---------------------------------------------------
form <- y ~ x * z + (x | species)
parsedForm <- glmercFormula(form, df, covList = covList, strList = list())


## ----, update parsed formula---------------------------------------------
covarSim <- c(0.5, -0.2, 0.5)
parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(covarSim))


## ----, update simulations------------------------------------------------
X <- model.matrix(nobars(form), df) # fixed effects design matrix
Z <- t(parsedForm$Lambdat %*% parsedForm$Zt) # random effects design
                                             # matrix with
                                             # phylogenetic
                                             # covariances
fixefSim <- rnorm(ncol(X)) # fixed effects
u <- rnorm(ncol(Z)) # whitened random effects
p <- plogis(as.numeric(X %*% fixefSim + Z %*% u)) # probability of observation
dl$y <- rbinom(nrow(df), 1, p) # presence-absence data
df <- as.data.frame(dl) # reconstruct the data frame with new
                        # structured response


## ----, plot lambda, fig.width=3, fig.height=3----------------------------
image(parsedForm$Lambdat)
image(crossprod(parsedForm$Lambdat))


## ----, plot Zt, fig.width=7, fig.height=2--------------------------------
image(parsedForm$Zt)


## ----, plot full cov, fig.width=5, fig.height=5--------------------------
image(fullCov <- t(parsedForm$Zt) %*% crossprod(parsedForm$Lambdat) %*% parsedForm$Zt)


## ----, plot sim data, fig.width=3, fig.height=3--------------------------
color2D.matplot(dl$y, xlab = "species", ylab = "sites", main = "abundance")


## ----, fit phylo model---------------------------------------------------
(mod <- glmerc(form, df, covList = covList))


## ----, compare fit and sim-----------------------------------------------
cbind(estimated = mod$opt$par, # estimated parameters
      true = c(covar = covarSim, fixef = fixefSim)) # true parameters


## ------------------------------------------------------------------------
data(fish)
data(limn)
Y <- as.matrix(fish)
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))
dl <- data.list(Y = t(Y), x = x,
                dimids = c("species", "sites"))
summary(dl)


