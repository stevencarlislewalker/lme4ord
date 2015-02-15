library(lme4ord)

n <- 10
m <- 30
dl <- dims_to_vars(data.list(Y = 1 * (rmat(n, m) > 0), x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)

modMat  <- list(model.matrix( ~ x, df),
                model.matrix( ~ 0 + z, df))
grpFac1 <- list(df$sites, df$sites)
grpFac2 <- list(df$species, df$species)
covMat1 <- lapply(rep(list(crossprod(rmat(n + 1, n))), 2), cov2cor)
covMat2 <- lapply(rep(list(crossprod(rmat(m + 1, m))), 2), cov2cor)

for(i in 1:3) dev.new()

ret <- mkTemplateReTrms(modMat, grpFac1, grpFac2, covMat1, covMat2)

X <- do.call(cbind, modMat)
Z <- t(ret$Lambdat %*% ret$Zt)
beta <- rnorm(ncol(X))
u <- rnorm(ncol(Z))
dl$Y <- rbinom(nrow(df), 1, plogis(as.numeric(X %*% beta + Z %*% u)))

ret <- within(ret, Lambdat@x[] <- mapToCovFact(c(1, -0.5, 1, 1)))

dev.set(2)
image(ret$Lambdat)
#image(ret$Lambdat[1:12, 1:12])
#image(ret$Lambdat[13:18, 13:18])
dev.set(3)
image(ret$Zt)
dev.set(4)
image(ret$Lambdat %*% ret$Zt)

parInds <- list(covar = 1:4, fixef = 5, loads = NULL)
initPars <- c(covar = c(1, 0, 1, 1), fixef = 0)
dfun <- mkGeneralGlmerDevfun(df$Y, model.matrix(Y ~ 1, df),
                             ret$Zt, ret$Lambdat,
                             rep(1, nrow(df)), rep(0, nrow(df)),
                             initPars, parInds,
                             ret$mapToCovFact, function(loads) NULL)
dfun(initPars)

opt <- optim(initPars, dfun, method = "L-BFGS-B",
             lower = c(0, -Inf, 0, 0, -Inf),
             control = list(trace = 4L))








cbind(.bdiag(ret$LambdatTheta)@x,
      .bdiag(ret$LambdatLind)@x,
      ret$Lambdat@x,
      ret$theta[ret$Lind])

plot(cbind(.bdiag(ret$LambdatBaseline)@x, ret$Lambdat@x))
plot(ret$Lambdat@x, ret$LambdatxBaseline * ret$theta[ret$Lind])
plot(ret$Lambdat@x, .bdiag(ret$LambdatBaseline)@x * ret$theta[ret$Lind])

as.vector(outer(outer(paste("species", 1:3),
                      paste("variable", 1:2), paste),
                      paste("species", 1:2), paste))

names(ret)
ret$Lambda <- t(ret$Lambdat)

image(ret$Zt)

image(crossprod(ret$Lambdat))
0.4 * 0.8 - 0.2 * 0.2

with(ret, cbind(Lambdat@x, Lind, theta[Lind], LambdatxBaseline))
cbind(ret$Lambdat@x, ret$theta[ret$Lind] * ret$LambdatxBaseline)
abline(0, 1)

ret2 <- mkTemplateReTrm(modMat[[2]], grpFac2[[2]], grpFac1[[2]], covMat2[[2]], covMat1[[2]])
ret2 <- within(ret2, Lambdat@x[] <- LambdatxBaseline * theta[Lind])
image(ret2$Lambdat)


ret1 <- mkTemplateReTrm(modMat[[1]], grpFac2[[1]], grpFac1[[1]], covMat2[[1]], covMat1[[1]])
ret1 <- within(ret1, Lambdat@x[] <- LambdatxBaseline * theta[Lind])
image(ret1$Lambdat)
