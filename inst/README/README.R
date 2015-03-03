
## ----, echo = FALSE------------------------------------------------------
library(knitr)
opts_chunk$set(fig.path = "inst/README/figure/")


## ----, message = FALSE---------------------------------------------------
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


## ----, eval = FALSE------------------------------------------------------
## glmerc(y ~ x * z + (x | species), data,
##        covList = list(species = Vphy),
##        family = binomial)


## ------------------------------------------------------------------------
set.seed(10)
n <- 100
m <- 10000
dl <- dims_to_vars(data.list(y = 1 * (matrix(rnorm(n * m), n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)


## ------------------------------------------------------------------------
phy <- rtree(n = m)
phy <- compute.brlen(phy, method = "Grafen", power = 0.1)


## ------------------------------------------------------------------------
Vphy <- stanCov(vcv(phy))
dimnames(Vphy) <- rep(list(1:m), 2)


## ----, fig.width = 4, fig.height = 5-------------------------------------
plot(phy)
image(as(Vphy, "sparseMatrix"))


## ------------------------------------------------------------------------
covList <- list(species = Vphy)


## ------------------------------------------------------------------------
form <- y ~ x*z + (x | species)
parsedForm <- glmercFormula(form, df, covList = covList)


## ------------------------------------------------------------------------
covarSim <- c(0.5, -0.2, 0.5)
parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(covarSim))


## ------------------------------------------------------------------------
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


## ----, fig.width=7, fig.height=2-----------------------------------------
image(parsedForm$Zt)


## ----, fig.width=5, fig.height=5-----------------------------------------
image(fullCov <- t(parsedForm$Zt) %*% crossprod(parsedForm$Lambdat) %*% parsedForm$Zt)


## ----, fig.width=3, fig.height=3-----------------------------------------
color2D.matplot(dl$y, xlab = "species", ylab = "sites", main = "abundance")


## ------------------------------------------------------------------------
system.time({
    (mod <- glmerc(form, df, covMat = covMat))
})


## ------------------------------------------------------------------------
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


