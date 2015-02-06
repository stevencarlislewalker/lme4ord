library(Matrix)
library(lme4ord)
library(lme4)
library(lme4pureR)
library(multitable)
library(pryr)

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

dl <- data.list(Y = t(Y), z = z, x = x, latent = 1:m,
                dimids = c("species", "sites"))
summary(dl)

m <- gblmer(Y ~ x*z + (1 | sites) + (1 | species), . ~ 0 + (0 + latent | sites),
            dl, binomial, 1, 1)
m

fitY <- matrix(getME(m, "mu"), 50, 30, byrow = TRUE)
boxplot(fitY ~ Y)

plot(fitY[,4], Y[,4])

plot(u, getME(m, "b")[1:n])
plot(v, m@loadings)

j <- 30
b <- getME(m, "b")[1:n]
bb <- seq(min(b), max(b), length = 100)
plot(b, Y[,j])
lines(bb, plogis(m@loadings[j]*bb))
