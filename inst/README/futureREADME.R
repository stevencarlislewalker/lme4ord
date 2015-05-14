library(lme4ord)
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


## ------------------------------------------------------------------------
set.seed(10)
n <- 10
m <- 30
dl <- dims_to_vars(data.list(y = 1 * (rmat(n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)


## ------------------------------------------------------------------------
phy <- rtree(n = m)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)


## ------------------------------------------------------------------------
Vphy <- stanCov(vcv(phy))
dimnames(Vphy) <- rep(list(1:m), 2)


covListNest <- list(`sites:species` =  diag(1, n, n) %x% Vphy)
form <- y ~ x * z + (1 | sites:species)
(modNest <- glmerc(form, df, binomial, covListNest))
#(modNest <- glmerc(form, df, binomial))



modNest$opt$par
sqrt(diag(solve(0.5 * lme4:::deriv12(modNest$dfun, modNest$opt$par)$Hessian)))

VarCorr(modNest)$`sites:species`

image(as(covListNest[[1]], "sparseMatrix"))

rho <- environment(modNest$dfun)
image(rho$pp$Lambdat)
image(environment(modNest$dfun)$pp$Zt)
plot(environment(modNest$dfun)$pp$u(1),
     as.numeric(getMEc(modNest, "X") %*% fixef(modNest)))

names(modNest$parsedForm)
modNest$parsedForm$Lind

getModMatAndGrpFac(findbars(form)[[1]], df)

