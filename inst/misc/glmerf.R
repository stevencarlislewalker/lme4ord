library(reo)
library(lme4ord)
data(fish)
data(limn)
Y <- as.matrix(fish)
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))

dl <- data.list(Y = Y, x = x,
                dimids = c("sites", "species"))
summary(dl)

df <- as.data.frame(dims_to_vars(dl))

U <- matrix(rnorm(n * 2), n, 2)
V <- matrix(rnorm(m * 2), m, 2)
E <- matrix(rnorm(n * m), n, m)
p <- plogis((U %*% t(V)) + E)
Y <- matrix(rbinom(n * m, 1, p), n, m)
dl <- variable(Y, c("sites", "species"))

mod <- glmerf(Y ~ 1 + (1 | sites) + (1 | species), dl, binomial,
              verbose = 3L,
              latentDims = 2L, loadingsDim = 2L, loadingPen = 5L,
              control = list(maxfun = 1500))
eigen(mod@optinfo$derivs$Hessian)$values
mod@optinfo$derivs$gradient

cbind(loadings = loadings(mod)[1:59],
      se = sqrt(diag(solve(0.5 * mod@optinfo$derivs$Hessian)))[1:59])

plot(V[,1], loadings(mod)[,1])
plot(V[,2], loadings(mod)[,2])
plot(loadings(mod))
ci <- c(-1, 1) * 1.96 * mean(sqrt(diag(solve(0.5 * mod@optinfo$derivs$Hessian)))[1:59])
abline(h = ci, v = ci)

image(tcrossprod(loadings(mod)) + diag(VarCorr(mod)$sites.1[1], m, m))

loadings(mod)



with(getME(mod, c("Lambdat", "Zt")), {
    image(t(Zt) %*% t(Lambdat) %*% Lambdat %*% Zt)
})

image(as(tcrossprod(loadings(mod)), "sparseMatrix"))
