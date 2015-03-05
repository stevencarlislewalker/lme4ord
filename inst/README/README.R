
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
n <- 95 # 100
m <- 475 # 500
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
    (mod <- glmerc(form, df, covMat = covMat,
                   optControl = list(iprint = 4L, maxfun = 500)))
})

centDfun <- function(par) mod$dfun(par) - mod$opt$fval
hh <- lme4:::deriv12(centDfun, mod$opt$par)
eigen(hh$Hessian)$values
1.96 * sqrt(diag(solve(0.5*hh$Hessian)))

library(numDeriv)
hhRich <- hessian(centDfun, mod$opt$par)
eigen(hh$Hessian)$values
cov2cor(solve(0.5 * hh$Hessian))
1.96 * sqrt(diag(solve(0.5*hhRich)))

centDfun

cbind(estimated = mod$opt$par, # estimated parameters
      true = c(covar = covarSim, fixef = fixefSim)) # true parameters



vc <- vcov(mod, FALSE)

1.96 * (sqrt(diag(vc))[1:3]/covar(mod))
log(covar(mod))

optPar <- mod$opt$par
pfun <- local({
    whichPar <- 1
    op <- optPar[-whichPar]
    function(par) {
        fn <- function(op) {
            parNow <- numeric(length(op) + 1)
            parNow[whichPar] <- par
            parNow[-whichPar] <- op
            centDfun(parNow)
        }
        opt <- bobyqa(op, fn, lower = mod$lower[-whichPar],
                      control = list(iprint = 4L, maxfun = 500))
        op <<- opt$par
        return(opt$fval)
    }
})

grd <- seq(mod$opt$par[whichPar] - 0.1, mod$opt$par[whichPar] + 0.1, length = 10)
pro <- sapply(grd, pfun)
plot(grd, pro, type = "o")

mod2 <- bobyqa(mod$opt$par, mod$dfun, lower = mod$lower,
               control = list(iprint = 4L))
mod2 <- bobyqa(mod2$par, mod$dfun, lower = mod$lower,
               control = list(iprint = 4L))

mod$opt$par
mod2$par

plot(interpSpline(grd, pro))
(cs <- c(0, qchisq(0.95, 1)))
abline(h = cs)
abline(v = mod$opt$par[whichPar])
abline(v = c(covarSim, fixefSim)[whichPar])
abline(v = mod$opt$par[whichPar] + c(-1, 1) * 1.96 * sqrt(diag(solve(0.5*hh$Hessian)))[whichPar])
xx <- seq(min(grd), max(grd), length = 100)
mQuad <- lm(pro ~ grd + I(grd^2))
lines(xx, predict(mQuad, newdata = data.frame(grd = xx)))

interpSpline(grd, pro)


interpSpline(grd, pro)$coef
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


