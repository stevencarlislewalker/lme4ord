Notes on `glmer` deviance functions
===================================


```r
library("lme4")
```

```
## Loading required package: Matrix
## Loading required package: methods
## Loading required package: Rcpp
```


```r
parsedForm <- glFormula(cbind(incidence, size - incidence) ~ period + (1 | herd),
                        family = binomial, data = cbpp)
```

`glmer` models are usually fitted in two steps.  The first step
involves no integration (`nAGQ = 0`), and the second involves
integration by either Laplace (`nAGQ = 1`) or adaptive Gaussian
quadrature (`nAGQ > 1`).  


```r
(devfun0 <- do.call(mkGlmerDevfun, parsedForm))
```

```
## function (theta) 
## {
##     resp$updateMu(lp0)
##     pp$setTheta(theta)
##     p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L), compDev, 
##         verbose)
##     resp$updateWts()
##     p
## }
## <environment: 0x7fbd1a031cd8>
```


```r
(devfun1 <- updateGlmerDevfun(devfun0, parsedForm$reTrms))
```

```
## function (pars) 
## {
##     resp$setOffset(baseOffset)
##     resp$updateMu(lp0)
##     pp$setTheta(as.double(pars[dpars]))
##     spars <- as.numeric(pars[-dpars])
##     offset <- if (length(spars) == 0) 
##         baseOffset
##     else baseOffset + pp$X %*% spars
##     resp$setOffset(offset)
##     p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, 
##         verbose)
##     resp$updateWts()
##     p
## }
## <environment: 0x7fbd1a031cd8>
```


```r
as.list(body(devfun1)[c(3, 4, 8, 9)])
```

```
## [[1]]
## resp$updateMu(lp0)
## 
## [[2]]
## pp$setTheta(as.double(pars[dpars]))
## 
## [[3]]
## p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose)
## 
## [[4]]
## resp$updateWts()
```

#### First stage deviance function

Note that there are several objects being used that do not appear to
be created anywhere (e.g. `lp0`).  Don't worry, they are actually in the
environment of `devfun0`, which we can list using the following
command.

```r
ls(rho <- environment(devfun0))
```

```
##  [1] "baseOffset"  "compDev"     "dpars"       "fac"         "GQmat"      
##  [6] "lower"       "lp0"         "nAGQ"        "pp"          "pwrssUpdate"
## [11] "resp"        "tolPwrss"    "verbose"
```
Therefore, all of these objects are available to `devfun0`.

We will now step through each line, explaining what it does.  The
first line,

```
## resp$updateMu(lp0)
```
updates the conditional mean of the response variable at the value of
the linear predictor, `lp0`.  Mathematically, the linear predictor is,
$$ \mathbf\eta = \mathbf o + \mathbf X\mathbf\beta + \mathbf Z
\mathbf\Lambda_\theta \mathbf u $$
The conditional mean is then,
$$ \mathbf\mu = \text{link}(\mathbf\eta) $$
Note that the condition mean depends on the covariance parameters
$\mathbf\theta$.  Therefore, this update of the conditional mean is
based on the old value of $\mathbf\theta$.  The new value of
$\mathbf\theta$ is updated on the next line:

```
## pp$setTheta(theta)
```
Next the penalized weighted residual sum of squares is calculated with
this new value of $\mathbf\theta$:

```
## p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L), compDev, verbose)
```
Finally the weights are updated (for safety?)

```
## resp$updateWts()
```


------------------------------------------------------------

#### Second stage deviance function




```r
ls(rho <- environment(devfun1))
```

```
##  [1] "baseOffset"  "compDev"     "dpars"       "fac"         "GQmat"      
##  [6] "lower"       "lp0"         "nAGQ"        "pp"          "pwrssUpdate"
## [11] "resp"        "tolPwrss"    "verbose"
```


```r
rho$baseOffset
```

```
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [36] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```



```
## resp$setOffset(baseOffset)
```


```
## resp$updateMu(lp0)
```


```
## pp$setTheta(as.double(pars[dpars]))
```


```
## spars <- as.numeric(pars[-dpars])
```


```
## offset <- if (length(spars) == 0) baseOffset else baseOffset + 
##     pp$X %*% spars
```


```
## resp$setOffset(offset)
```



