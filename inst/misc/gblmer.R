library(Matrix)
library(lme4ord)
library(lme4)
library(lme4pureR)
library(multitable)
library(pryr)
library(reo)

n <- 50
m <- 30
x <- rnorm(n)
z <- rnorm(m)
B <- rmat(2, 2)
u <- rnorm(n)
v <- rnorm(m)

eta <-
    (cbind(1, x) %*% B %*% rbind(1, z)) +
    (cbind(1, x) %*%        rmat(2, m)) +
    ( rmat(n, 2) %*%       rbind(1, z)) +
    (   cbind(u) %*%          rbind(v))

Y <- rmat(n, m,
          f(rbinom(n, size = 1, prob = plogis(eta))),
          eta = eta)


data(fish)
data(limn)
Y <- as.matrix(fish)
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))

dl <- data.list(Y = Y, x = x,
                dimids = c("sites", "species"))
summary(dl)

mod00 <- gblmer(Y ~ 1 + (1 | species), dl, binomial,
                loadingsDim = 0)

mod01 <- gblmer(Y ~ x + (x | species), dl, binomial,
                loadingsDim = 0)

mod0 <- gblmer(Y ~ 1 + (1 | species), dl, binomial,
               loadingsDim = 2,
               latentDim = 2,
               loadingPen = 1L,
               verbose = 2L)

mod1 <- gblmer(Y ~ x + (x | species), dl, binomial,
               loadingsDim = 2,
               latentDim = 2,
               loadingPen = 1L,
               verbose = 2L)

mod00
mod01
mod0
mod1
summary(mod0)
summary(mod1)

mod <- mod0
ranef(mod)$species
ranef(mod)$sites
(loadMod <- loadings(mod))
latentCov <- Matrix(loadMod %*% t(loadMod)) + diag(VarCorr(mod)$species[1], m, m)
image(cov2cor(latentCov))
image(vcov(mod)@factors$correlation)

mods <- list(mod00, mod01, mod0, mod1)
fits <- lapply(mods, function(mm) matrix(getME(mm, "mu"), n, m))
par(mfrow = c(2, 2))
trash <- lapply(fits, function(ff) boxplot(ff ~ Y, las = 1, horizontal = TRUE))

k <- "latent1"
mod <- mod0
b <- ranef(mod)$sites[, k]
bb <- seq(min(b), max(b), length = 100)
for(j in 1:m) {
    plot(ranef(mod)$sites[ , k], Y[, j],
         main = colnames(Y)[j])
    lines(bb, plogis(loadings(mod)[j, k]*bb))
    readline("press any key")
}

j <- 3
k <- "latent2"
mod <- mod1
b <- ranef(mod)$sites[, k]
bb <- seq(min(b), max(b), length = 100)
plot(b, Y[,j])
lines(bb, plogis(loadings(mod)[j, k]*bb))

par(mfrow = c(1, 2))
for(mod in c(mod0, mod1)) {
    vm <- varimax(loadings(mod))
    biplot(as.matrix(ranef(mod)$sites) %*% vm$rotmat,
           vm$loadings, cex = 0.5)
}


plot(loadings(mod0), loadings(mod1))
