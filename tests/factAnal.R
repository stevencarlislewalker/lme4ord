library(lme4ord)
library(multitable)
library(reo)

data(fish)
data(limn)
dataList <- dims_to_vars(data.list(respVar = as.matrix(fish),
                                   pH = limn$pH,
                                   dimids = c("lakes", "species")))
dataFrame <- as.data.frame(aperm(dataList, c(2, 1)))

gm1 <- strucGlmer(respVar ~ 1 + 
                  (1 | lakes) +
                  (1 | species) +
                  factAnal(0 + lakes | species, nAxes = 2),
                  family = binomial(), data = dataFrame,
                  devfunOnly = FALSE,
                  optMaxit = 20000, optVerb = 0L,
                  penLoads = mkPenLpNorm(p = 2))

inds <- reIndsPerTrm(gm1)
with(gm1$parsedForm, image(tcrossprod(Lambdat * Zt)[inds[[2]], inds[[3]]]))
with(gm1$parsedForm, image(tcrossprod(Lambdat * Zt)))
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut))
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[1]], inds[[1]]])
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[1]], inds[[2]]])
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[1]], inds[[3]]])
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[2]], inds[[2]]])
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[2]], inds[[3]]])
image(tcrossprod(gm1$parsedForm$devfunEnv$pp$Ut)[inds[[3]], inds[[3]]])

gm2 <- strucGlmer(respVar ~ 1 + 
                  (1 | lakes) +
                  (1 | species) +
                  factAnal(0 + lakes | species, nAxes = 2, seed = 1),
                  family = binomial(), data = dataFrame,
                  devfunOnly = FALSE,
                  optMaxit = 20000, optVerb = 0L,
                  penLoads = mkPenLpNorm())

loadMat1 <- as.matrix(repSparseGenFullTri(30, 2, loads(gm1)))
loadMat2 <- as.matrix(repSparseGenFullTri(30, 2, loads(gm2)))

## correct sign flips
Q <- orthProcrustesRotMat(loadMat1, loadMat2)
loadMat1 <- loadMat1 %*% Q
# check that rotation is _only_ corrects sign flips, ie that:
                                        # abs(diagonal) is one
all.equal(abs(diag(Q)), rep(1, nrow(Q)), tolerance = 1e-5)
                                        # off-diagonal is zero
all.equal(Q[c(2, 3)],   rep(0, nrow(Q)), tolerance = 1e-5)


plot(loadMat1[,1], loadMat2[,1])
abline(a = 0, b = 1)

plot(loadMat1[,2], loadMat2[,2])
abline(a = 0, b = 1)

apply(loadMat1, 2, mean)
apply(loadMat2, 2, mean)

apply(loadMat1, 2, sd)
apply(loadMat2, 2, sd)

plot(loads(gm1), gm1$opt$par[-(1:3)])
plot(loads(gm2), gm2$opt$par[-(1:3)])

axes1 <- matrix(ranef(gm1)[[3]], 52, 2) %*% Q
axes2 <- matrix(ranef(gm2)[[3]], 52, 2)

par(mfrow = c(1, 2))
biplot(axes1[,2:1], loadMat1[,2:1])
biplot(axes2[,2:1], loadMat2[,2:1])

gm1$parsedForm$Zt
gm2$parsedForm$Zt

plot(ranef(gm1)[[1]], ranef(gm2)[[1]])
abline(a = 0, b = 1)
plot(ranef(gm1)[[2]], ranef(gm2)[[2]])
abline(a = 0, b = 1)

plot(fitted(gm1), fitted(gm2))
abline(a = 0, b = 1)

## non-standard models
gm0 <- strucGlmer(respVar ~ (1 | lakes) + (1 | species) +
                  factAnal(0 + lakes | species, nAxes = 1),
                  family = binomial(), data = dataFrame,
                  devfunOnly = FALSE,
                  optMaxit = 20000, optVerb = 0L,
                  penLoads = mkPenLpNorm())

gm0pH <- strucGlmer(respVar ~ (1 | lakes) + (1 | species) +
                    factAnal(pH | species, nAxes = 1),
                    family = binomial(), data = dataFrame,
                    devfunOnly = FALSE,
                    optMaxit = 20000, optVerb = 0L,
                    penLoads = mkPenLpNorm())
