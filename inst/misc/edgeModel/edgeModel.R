library(lme4ord)
library(lme4)
library(ape)

td <- simTestPhyloDat(9, n = 10, m = 100,
                      form = y ~ 1 + (1 | species),
                      power = 0.1)
getReTrm.edgeStruct <- function(){}


form <- y ~ x * z + edgeStruct(1 | species, ph = ph) + 1 + (x | sites)
splitForm(form)$reTrmAddArgs[[1]]





generalParseFormula(form, as.data.frame(td$dl))

edgeInd <- edgeTipIndicator(td$ph)
dummy <- as.data.frame(t(edgeInd))
td$dl <- td$dl + variableGroup(dummy, "species")
edgeNms <- names(dummy)
df <- as.data.frame(td$dl)

image(as(tcrossprod(edgeInd), "sparseMatrix"))

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

opt <- optim(c(1, 0), dfun, lower = c(0, -Inf), method = "L-BFGS-B")
dfun(opt$par)
opt$par

rho <- environment(dfun)
with(rho$pp, image(crossprod(Lambdat %*% Zt)[1:50, 1:50]))

plot(td$ph)
edgelabels(round(rho$pp$b(1), 2), cex = 0.6)
