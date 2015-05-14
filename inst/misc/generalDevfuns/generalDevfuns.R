
## ----prelim, message = FALSE, results = 'hide', echo = FALSE, cache = FALSE----
library("lme4ord")
library("lme4")
library("multitable")
opts_chunk$set(cache = TRUE)


## ----simLmer-------------------------------------------------------------
form <- y ~ 1 + (0 + x:g1 | g2)
dl <- data.list(y = matrix(rnorm(10 * 1), 10, 1), x = rnorm(10), 
                g1 = letters[1:10], g2 = LETTERS[1],
                dimids = c("g1", "g2"))
df <- as.data.frame(dl)
m <- lmer(form, df, control = lmerControl(
                        check.nobs.vs.nRE = "ignore",
                        check.nlev.gtr.1 = "ignore"))


## ----simulate params-----------------------------------------------------
n <- 60 # samples 
p <- 2  # fixed effects
q <- 30 # random effects


## ----cov fact sims-------------------------------------------------------
covTemplate <- as(chol(cov(matrix(rnorm((q+1)*q), q + 1, q))),
                  "sparseMatrix")
Lambdat <- covTemplate


## ----viewLambdat, fig.width = 4, fig.height = 4--------------------------
image(Lambdat)


## ----mod mat sims--------------------------------------------------------
Zt <- sparseMatrix(i = rep(1:q, 6),
                   j = as.vector(outer(rep(1:10, 3), seq(0, 50, 10), "+")),
                   x = rep(rnorm(10), 18))
X <- matrix(rnorm(n * p), n, p)


## ----viewZt, fig.width = 5, fig.height = 3-------------------------------
image(Zt)


## ----resp vect sims------------------------------------------------------
eta <- as.numeric(X %*% c(1, 1) + t(Zt) %*% t(Lambdat) %*% rnorm(q))
y <- rbinom(n, 1, plogis(eta))
weights <- rep(1, n); offset <- rep(0, n)


## ----initial-------------------------------------------------------------
initPars <- c(covar = 1,
              fixef = c(0, 0),
              loads = rnorm(10))


## ----par indices---------------------------------------------------------
parInds <- list(covar = 1,
                fixef = 2:3,
                loads = 4:13)


## ----mappings------------------------------------------------------------
mapToCovFact <- function(covar) covar * covTemplate@x
mapToModMat <- function(loads) rep(loads, 18)


## ----illustrate mappings-------------------------------------------------
Lambdat@x <- mapToCovFact(initPars[parInds$covar])
Zt@x <- mapToModMat(initPars[parInds$loads])


## ----make the deviance function------------------------------------------
devfun <- mkGeneralGlmerDevfun(y, X, Zt, Lambdat,
                               weights, offset,
                               initPars, parInds,
                               mapToCovFact, mapToModMat)                               


## ----use the dev funct---------------------------------------------------
devfun(initPars)


## ----optimize the deviance function--------------------------------------
opt <- minqa:::bobyqa(initPars, devfun)


## ----look at the optimal solution----------------------------------------
setNames(opt$par, names(initPars))


## ----more packages, message = FALSE--------------------------------------
library(reo)
library(multitable)


## ----fish data-----------------------------------------------------------
data(fish)
data(limn)
Y <- as.matrix(fish)
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))
dl <- data.list(Y = Y, x = x,
                dimids = c("sites", "species"))
dl <- dims_to_vars(dl)
summary(dl)


## ----long----------------------------------------------------------------
head(df <- as.data.frame(dl))


## ----fish resp-----------------------------------------------------------
y <- df$Y
weights <- rep(1, length(y)); offset <- rep(0, length(y))


## ----fish X--------------------------------------------------------------
X <- model.matrix(Y ~ x, df)[,]


## ----fish Zt-------------------------------------------------------------
Jspecies <- t(as(as.factor(df$species), Class = "sparseMatrix"))
Zt <- KhatriRao(t(Jspecies), t(X))


## ----loadings matrix-----------------------------------------------------
nFreeLoadings <- m
U <- matrix(1:nFreeLoadings, nrow = m, ncol = 1)
latentVarNames <- "latent"
U <- setNames(as.data.frame(U), latentVarNames)
latentData <- data.list(U, drop = FALSE, dimids = "species")
df <- as.data.frame(dl + latentData)
Jsites <- with(df, t(as(as.factor(sites), Class = "sparseMatrix")))
Zt <- rBind(KhatriRao(t(Jsites), t(model.matrix(Y ~ 0 + latent, df))), Zt)


## ----cov fac fish--------------------------------------------------------
templateFact <- sparseMatrix(i = 1:n, j = 1:n, x = rep(1, n))
ii <- rep(1:2, 1:2); jj <- sequence(1:2)
templateRanef <- sparseMatrix(i = ii, j = jj, x = 1 * (ii == jj))
Lambdat <- .bdiag(c(list(templateFact), 
                    rep(list(templateRanef), m)))


## ----loadings mapping fish-----------------------------------------------
mapToModMat <- local({
    Ztx <- Zt@x
    Zwhich <- Zt@i %in% (seq_len(n) - 1)
    Zind <- Zt@x[Zwhich]
    loadInd <- 1:nFreeLoadings
    function(loads) {
        Ztx[Zwhich] <- loads[Zind]
        return(Ztx)
    }
})


## ----lambda mapping fish-------------------------------------------------
mapToCovFact <- local({
    Lambdatx <- Lambdat@x
    LambdaWhich <- (n+1):length(Lambdatx)
    LindTemplate <- ii + 2 * (jj - 1) - choose(jj, 2)
    Lind <- rep(LindTemplate, m)
    function(covar) {
        Lambdatx[LambdaWhich] <- covar[Lind]
        return(Lambdatx)
    }
})


## ----pars fish-----------------------------------------------------------
initPars <- c(covar = c(1, 0, 1),
              fixef = c(0, 0),
              loads = rep(0, m))
parInds <- list(covar = 1:3,
                fixef = 4:5,
                loads = (1:m)+5)


## ----make fish devfun, warning = FALSE-----------------------------------
devfun <- mkGeneralGlmerDevfun(y, X, Zt, Lambdat,
                               weights, offset, 
                               initPars, parInds,
                               mapToCovFact, mapToModMat)


## ----opt fish, eval = TRUE-----------------------------------------------
opt <- minqa:::bobyqa(initPars, devfun)
setNames(opt$par, names(initPars))


## ----ordination diagram, fig.width = 7, fig.height = 7-------------------
siteScores <- cbind(pH = x,
                    latent = environment(devfun)$pp$b(1)[1:n])
speciesScores <- cbind(pH = environment(devfun)$pp$b(1)[-(1:(n+m))],  
                       latent = opt$par[parInds$loads])
rownames(siteScores) <- rownames(Y)
rownames(speciesScores) <- colnames(Y)
biplot(siteScores, speciesScores, cex = 0.5)


