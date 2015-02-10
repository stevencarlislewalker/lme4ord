lme4ord (l-m-e-ford)
====================



Mixed-effects models for community ecologists.  See the currently
evolving [mission statement](https://github.com/stevencarlislewalker/lme4ord/issues/1).

This package is not at all stable.

#### Newer function: `mkGeneralGlmerDevfun`

$$ \sin(x) = \sqrt{1 - 2x} $$

#### New function: `gblmer`


```r
library(Matrix)
```

```
## Loading required package: methods
```

```r
library(lme4ord)
```

```
## Loading required package: lme4
## Loading required package: Rcpp
## 
## Attaching package: 'lme4ord'
## 
## The following object is masked from 'package:stats':
## 
##     loadings
```

```r
library(lme4)
library(lme4pureR)
library(multitable)
library(pryr)
```

```
## Warning: package 'pryr' was built under R version 3.1.2
```

```r
library(reo)
```

```
## Loading required package: MASS
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
## Loading required package: ellipse
```

```r
data(fish)
data(limn)
Y <- as.matrix(fish)
## Y <- Y[, colSums(Y) > 1]
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))
dl <- data.list(Y = t(Y), x = x,
                dimids = c("species", "sites"))
summary(dl)
```

```
##            Y     x
## species TRUE FALSE
## sites   TRUE  TRUE
```

```r
mod <- gblmer(Y ~ 1 + (1 | species), . ~ 0 + (0 + latent | sites),
              dl, binomial, 1, 1, 2)
```

```
## Error in gblmer(Y ~ 1 + (1 | species), . ~ 0 + (0 + latent | sites), dl, : data must be a data list
```

```r
mod
```

```
## Error in eval(expr, envir, enclos): object 'mod' not found
```

```r
ranef(mod)$species
```

```
## Error in ranef(mod): object 'mod' not found
```

```r
ranef(mod)$sites
```

```
## Error in ranef(mod): object 'mod' not found
```

```r
(loadMod <- loadings(mod))
```

```
## Error in loadings(mod): object 'mod' not found
```

```r
latentCov <- Matrix(loadMod %*% t(loadMod)) + diag(VarCorr(mod)$species[1], m, m)
```

```
## Error in is(data, "Matrix"): object 'loadMod' not found
```

```r
image(cov2cor(latentCov))
```

```
## Error in image(cov2cor(latentCov)): error in evaluating the argument 'x' in selecting a method for function 'image': Error in cov2cor(latentCov) : 
##   error in evaluating the argument 'V' in selecting a method for function 'cov2cor': Error: object 'latentCov' not found
```

```r
fitY <- matrix(getME(mod, "mu"), n, m, byrow = TRUE)
```

```
## Error in is(object, "merMod"): object 'mod' not found
```

```r
boxplot(fitY ~ Y, las = 1, horizontal = TRUE)
```

```
## Error in eval(expr, envir, enclos): object 'fitY' not found
```

#### Short demo

We need these packages (and their dependencies),

```r
library("lme4")
library("lme4ord")
library("Matrix")
library("reo") ## install.packages("reo", repos="http://R-Forge.R-project.org")
```
Prepare some data,

```r
Y <- Yp <-  as.matrix(fish)
Y <- Y[order(rowSums(Y)), ]
Yp <- Yp[order(rowSums(Yp)), ]
```
FIXME:  use more standard data set in more standard package.

Construct deviance functions for one and two axis ordination models,

```r
dfun1 <- logisticPcaDevfun(Yp, 1)
```

```
## Note: method with signature 'Matrix#diagonalMatrix' chosen for function 'kronecker',
##  target signature 'dgCMatrix#ddiMatrix'.
##  "sparseMatrix#ANY" would also be valid
## Note: method with signature 'dsparseMatrix#dsparseMatrix' chosen for function 'kronecker',
##  target signature 'dgCMatrix#dtTMatrix'.
##  "sparseMatrix#TsparseMatrix" would also be valid
```

```r
dfun2 <- logisticPcaDevfun(Yp, 2)
```
Get starting values for the optimization of these deviance functions,

```r
pars1 <- unlist(as.list(environment(dfun1))[c("theta", "phi")])[-1]
pars2 <- unlist(as.list(environment(dfun2))[c("theta", "phi")])[-1]
```
Optimize these deviance functions,

```r
opt1 <- optim(pars1, dfun1, method = "BFGS",
              control = list(maxit = 500, trace = TRUE))
```

```
## initial  value 1261.432956 
## iter  10 value 1169.299947
## iter  20 value 1168.688341
## iter  30 value 1168.118807
## iter  40 value 1167.438431
## final  value 1167.329338 
## converged
```

```r
opt2 <- optim(pars2, dfun2, method = "BFGS",
              control = list(maxit = 500, trace = TRUE))
```

```
## initial  value 1189.855260 
## iter  10 value 1105.507937
## iter  20 value 1102.367337
## iter  30 value 1101.121497
## iter  40 value 1099.596968
## iter  50 value 1098.988855
## iter  60 value 1098.701727
## iter  70 value 1098.460386
## iter  80 value 1098.395093
## iter  90 value 1098.384686
## iter  90 value 1098.384685
## final  value 1098.383082 
## converged
```
Both models seem to both converge,

```r
opt1$convergence
```

```
## [1] 0
```

```r
opt2$convergence
```

```
## [1] 0
```

FIXME:  However for the two-axis model, while it gets close quickly, takes hundreds of iterations zeroing in on a solution.  Is this a quasi-convex problem?  That is, is it just skating around on a fairly flat part of the deviance function?  Perhaps we could get a speed up with some kind of penalty?


We make easier to understand objects from the results,

```r
mod1 <- mkMod(environment(dfun1), opt1)
mod2 <- mkMod(environment(dfun2), opt2)
```

Let's plot some results from the two-axis model.  First we plot a series of image plots of observed and fitted site-by-species matrices.  These plots provide a decomposition of the sources of variation in the observed sites by species matrix (FIXME: add residual plot too).

```r
plotimage <- function(mat, ...)
    image(1:nrow(mat), 1:ncol(mat), mat, las = 1,
          zlim = c(0, 1),
          col = grey(seq(1, 0, length = 100)),
          ...)
par(mfrow = c(1, 6))
plotimage(Yp, main = "data")
plotimage(plogis(mod2$fit),
          main = "fitted values")
plotimage(plogis(mod2$fitInter),
          main = "intercept")
plotimage(plogis(mod2$fitAxes),
          main = "site-species interactions")
plotimage(plogis(mod2$fitRow),
          main = "main site effect")
plotimage(plogis(mod2$fitCol),
          main = "main species effect")
```

![plot of chunk unnamed-chunk-11](inst/README/figure/unnamed-chunk-11-1.png) 

Now we make a logit-scale biplot (with only a few species to reduce clutter),

```r
par(mfrow = c(1, 1))
rowKeep <- apply(abs(mod2$rowScores) > 0, 1, any)
colKeep <- apply(abs(mod2$colScores) > 0.3, 1, any)
biplot(mod2$rowScores[rowKeep,c(1, 2)], mod2$colScores[colKeep,c(1, 2)],
       xlabs = (1:52)[rowKeep], ylabs = colnames(Yp)[colKeep],
       xlab = "Axis I", ylab = "Axis II")
```

![plot of chunk unnamed-chunk-12](inst/README/figure/unnamed-chunk-12-1.png) 

Note that the two kinds of bass (smallmouth, SB, and largemouth, LB) are orthogonal, indicating that they are relatively uncorrelated.  On the other hand, northern redbelly dace, NRD, is negatively correlated with largemouth.

We can also plot the covariance matrix among species of the latent variables,

```r
image(cov2cor(mod2$typeCors))
```

![plot of chunk unnamed-chunk-13](inst/README/figure/unnamed-chunk-13-1.png) 

#### TODO

Lots!  Most important things:

1. write up math
2. user interface
3. allow arbitrary family
4. find faster parameterizations
