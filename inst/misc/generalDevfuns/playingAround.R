library(lme4ord)
library(multitable)
library(ape)

# Generate simulated data for nspp species and nsite sites
nspp <- 15
nsite <- 10

# residual variance (set to zero for binary data)
sd.resid <- 0

# fixed effects
beta0 <- 0
beta1 <- 0

# magnitude of random effects
sd.B0 <- 1
sd.B1 <- 1

# whether or not to include phylogenetic signal in B0 and B1
signal.B0 <- TRUE
signal.B1 <- TRUE

# simulate a phylogenetic tree
phy <- rtree(n = nspp)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

# standardize the phylogenetic covariance matrix to have determinant 1
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

# Generate environmental site variable
X <- matrix(1:nsite, nrow = 1, ncol = nsite)
X <- (X - mean(X))/sd(X)

# Perform a Cholesky decomposition of Vphy. This is used to
# generate phylogenetic signal: a vector of independent normal random
# variables, when multiplied by the transpose of the Cholesky
# deposition of Vphy will have covariance matrix equal to Vphy.

iD <- t(chol(Vphy))

# Set up species-specific regression coefficients as random effects
if (signal.B0 == TRUE) {
		b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
} else {
		b0 <- beta0 + rnorm(nspp, sd = sd.B0)
}
if (signal.B1 == TRUE) {
		b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
} else {
		b1 <- beta1 + rnorm(nspp, sd = sd.B1)
}

# Simulate species abundances among sites to give matrix Y that
# contains species in rows and sites in columns
y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
e <- rnorm(nspp * nsite, sd = sd.resid)
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)

Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
Y <- matrix(Y, nrow = nspp, ncol = nsite)

# name the simulated species 1:nspp and sites 1:nsites
rownames(Y) <- 1:nspp
colnames(Y) <- 1:nsite

par(mfrow = c(3, 1), las = 1, mar = c(2, 4, 2, 2) - 0.1)
matplot(t(X), type = "l", ylab = "X", main = "X among sites")
hist(b0, xlab = "b0", main = "b0 among species")
hist(b1, xlab = "b1", main = "b1 among species")

par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
if(require(plotrix))
    color2D.matplot(Y, ylab = "species", xlab = "sites", main = "abundance")

# Transform data matrices into "long" form, and generate a data frame
YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)

XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
nspp * nsite, ncol = 1)

site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
1)), nrow = nspp * nsite, ncol = 1)
sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
nrow = nspp * nsite, ncol = 1)

dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))

# Format input and perform communityPGLMM()
# set up random effects

# random intercept with species independent
re.1 <- list(1, sp = dat$sp, covar = diag(nspp))

# random intercept with species showing phylogenetic covariances
re.2 <- list(1, sp = dat$sp, covar = Vphy)

# random slope with species independent
re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))

# random slope with species showing phylogenetic covariances
re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)

# random effect for site
re.site <- list(1, site = dat$site, covar = diag(nsite))

dl <- dims_to_vars(data.list(Y = t(Y), X = as.numeric(X), const = rep("const", nrow(Y)),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
df$sites <- factor(df$sites, levels = unique(as.character(df$sites)))
df$species <- factor(df$species, levels = unique(as.character(df$species)))

y <- df$Y
weights <- rep(1, length(y)); offset <- rep(0, length(y))

M <- model.matrix(Y ~ X, df)

ordV <- order(rownames(Vphy))
VphyOrd <- Vphy[ordV, ordV]
dimnames(VphyOrd) <- rep(list(levels(df$species)), 2)

modMat <- c(rep(list(model.matrix(Y ~ 1, df)), 2),
            rep(list(model.matrix(Y ~ 0 + X, df)), 2),
            list(model.matrix(Y ~ 1, df)))
grpFac1 <- c(rep(list(df$species), 4), list(df$sites))
grpFac2 <- rep(list(df$const), 5)
covMat1 <- list(diag(nspp), VphyOrd, diag(nspp), VphyOrd, diag(nsite))
covMat2 <- rep(list(matrix(1, 1, 1)), 5)

ret <- mkTemplateReTrms(modMat, grpFac1, grpFac2, covMat1, covMat2)
image(ret$Zt)
image(ret$Lambdat)

initPars <- c(covar = rep(1, 5), fixef = rep(0, 2))
parInds <- list(covar = 1:5, fixef = 6:7, loads = NULL)
dfun <- mkGeneralGlmerDevfun(y, M, ret$Zt, ret$Lambdat,
                             rep(1, nspp * nsite), rep(0, nspp * nsite),
                             initPars, parInds,
                             ret$mapToCovFact, function(loads) NULL)
dfun(initPars)
opt <- optim(initPars, dfun, method = "L-BFGS-B",
             lower = c(rep(0, 5), rep(-Inf, 2)),
             control = list(trace = 0L))
dfun(opt$par)
image(environment(dfun)$pp$Lambdat)

## Jspecies <- as(df$species, "sparseMatrix")
## Jsites <- as(df$sites, "sparseMatrix")
## image(Jsites)
## image(Jspecies)
## image(KhatriRao(KhatriRao(Jsites, Jspecies), t(df$X)))
## image(KhatriRao(KhatriRao(Jspecies, Jsites), t(df$X)))
## image(KhatriRao(Jspecies, t(df$X)))
## image(Jspecies)
## image(Jsites)

# image(mkTemplateTermLambdat(re.4$covar, 10))
                            
## Ztlist <- list(KhatriRao(Jspecies, t(M[,1])),
##                KhatriRao(Jspecies, t(M[,1])),
##                KhatriRao(Jspecies, t(M[,2])),
##                KhatriRao(Jspecies, t(M[,2])))
## Zt <- do.call(rBind, Ztlist)

## Lambdatlist <- list(as(t(chol(re.1$covar)), Class="sparseMatrix"),
##                     as(t(chol(re.2$covar)), Class="sparseMatrix"),
##                     as(t(chol(re.3$covar)), Class="sparseMatrix"),
##                     as(t(chol(re.4$covar)), Class="sparseMatrix"))
## Lambdat <- .bdiag(Lambdatlist)

## Zt <- Zt[46:60, ]
## Lambdat <- Lambdat[46:60, 46:60]

## image(as(kronecker(matrix(1, 10, 15), t(chol(re.2$covar))), Class = "sparseMatrix"))

## image(Lambdat)
## image(Zt)
## image(Lambdat %*% Zt %*% t(Zt) %*% t(Lambdat))
## image(t(Zt) %*% t(Lambdat) %*% Lambdat %*% Zt)

## mapToModMat <- local({
##     Ztx <- Zt@x
##     function(loads) Ztx
## })

## mapToCovFact <- local({
    
## })
