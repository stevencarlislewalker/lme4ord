library("lme4ord")
library("multitable")
library("ape")
library("minqa")
library("tony")

## --------------------------------------------------------------------
## this example leads to lots of fits on the boundary.  but tony's
## code and lme4 give _similar_ answers.  we shouldn't expect _exact_
## answers because tony's code uses PQL and lme4ord uses the Laplace
## approximation
## --------------------------------------------------------------------

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
y <- matrix(outer(b0, array(1, dim = c(1, nsite))),
            nrow = nspp,
            ncol = nsite) +
    matrix(outer(b1, X),
           nrow = nspp,
           ncol = nsite)
e <- rnorm(nspp * nsite, sd = sd.resid)
y <- y + matrix(e, nrow = nspp, ncol = nsite)
y <- matrix(y, nrow = nspp * nsite, ncol = 1)
Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
Y <- matrix(Y, nrow = nspp, ncol = nsite)
# name the simulated species 1:nspp and sites 1:nsites
rownames(Y) <- 1:nspp
colnames(Y) <- 1:nsite
## par(mfrow = c(3, 1), las = 1, mar = c(2, 4, 2, 2) - 0.1)
## matplot(t(X), type = "l", ylab = "X", main = "X among sites")
## hist(b0, xlab = "b0", main = "b0 among species")
## hist(b1, xlab = "b1", main = "b1 among species")
## par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
## if(require(plotrix))
##     color2D.matplot(Y, ylab = "species", xlab = "sites", main = "abundance")
# Transform data matrices into "long" form, and generate a data frame
YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
nspp * nsite, ncol = 1)
site <- matrix(kronecker(1:nsite,
                         matrix(1, nrow = nspp, ncol = 1)),
               nrow = nspp * nsite, ncol = 1)
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
## re.site <- list(1, site = dat$site, covar = diag(nsite))
dl <- dims_to_vars(data.list(Y = t(Y), X = as.numeric(X),
                             dimids = c("sites", "species")))
dl <- dl + variable(dl$species, "species", "speciesInd")
df <- as.data.frame(dl)
df$sites <- factor(df$sites, levels = unique(as.character(df$sites)))
df$species <- factor(df$species, levels = unique(as.character(df$species)))
df$speciesInd <- factor(df$speciesInd, levels = unique(as.character(df$speciesInd)))

form <- Y ~ X +
    (1 | speciesInd) +
    (1 | species) +
    (0 + X | speciesInd) +
    (0 + X | species)
covList <- list(species = Vphy,
                speciesInd = diag(nspp))

modSteve <- glmerc(form, df, binomial, covList)


modTony <- communityPGLMM(Y ~ X, data = dat, family = "binomial",
                          sp = dat$sp, site = dat$site,
                          random.effects = list(re.1, re.2, re.3, re.4),
                          REML = FALSE, verbose = TRUE)

(covarComp <- cbind(tony = round(modTony$s2r, 3),
                    steve = round(modSteve$opt$par[modSteve$parInds$covar]^2, 3),
                    true = c(sd.B0, sd.B0, sd.B1, sd.B1)))
(fixefComp <- cbind(tony = as.vector(modTony$B),
                    steve = modSteve$opt$par[modSteve$parInds$fixef],
                    true = c(beta0, beta1)))
