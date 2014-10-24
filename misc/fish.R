library(lme4ord)
library(Matrix)
library(reo)

Y <- Yp <-  as.matrix(fish)
Y <- Y[order(rowSums(Y)), ]
Yp <- Yp[order(rowSums(Yp)), ]
dfun <- logisticPcaDevfun(Yp, 1, Upenalty = 0, thetaPenalty = 0)

pars <- unlist(as.list(environment(dfun))[c("theta", "phi")])[-1]
(opt <- optim(pars, dfun, method = "BFGS",
              control = list(maxit = 500, trace = TRUE)))


dfun(opt$par)
rho <- environment(dfun)
mod <- mkMod(environment(dfun), opt)

plotimage <- function(mat, ...)
    image(1:nrow(mat), 1:ncol(mat), mat, las = 1,
          zlim = c(0, 1),
          col = grey(seq(1, 0, length = 100)),
          ...)
par(mfcol = c(2, 5), mar = c(3, 3, 1, 1))
plotimage(Yp)
plotimage(plogis(mod$fit))
plotimage(Yp)
plotimage(plogis(mod$fitInter))
plotimage(Yp)
plotimage(plogis(mod$fitAxes))
plotimage(Yp)
plotimage(plogis(mod$fitRow))
plotimage(Yp)
plotimage(plogis(mod$fitCol))
par(mfrow = c(1, 1))
rowKeep <- apply(abs(mod$rowScores) > 0, 1, any)
colKeep <- apply(abs(mod$colScores) > 0.3, 1, any)
biplot(mod$rowScores[rowKeep,c(1, 2)], mod$colScores[colKeep,c(1, 2)],
       xlabs = (1:52)[rowKeep], ylabs = colnames(Yp)[colKeep],
       xlab = "Axis I", ylab = "Axis II")
biplot(mod$rowScores[rowKeep,c(1, 3)], mod$colScores[colKeep,c(1, 3)],
       xlabs = (1:52)[rowKeep], ylabs = colnames(Yp)[colKeep],
       xlab = "Axis I", ylab = "Axis III")

image(cov2cor(mod$typeCors))
