library("lme4ord")
library("minqa")
library("plotrix")
library("ape")

                                        # initial simulation of one
                                        # binary sites by species
                                        # matrix, with environmental
                                        # variable, x, and functional
                                        # trait, z (will change y
                                        # below to have structure)
set.seed(10)
n <- 10
m <- 30
dl <- dims_to_vars(data.list(y = 1 * (rmat(n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)

                                        # make up some silly
                                        # phylogeny
phy <- rtree(n = m)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)

                                        # estimate a phylogenetic
                                        # covariance matrix,
                                        # standardized to unit
                                        # determinant
Vphy <- stanCov(vcv(phy))
dimnames(Vphy) <- rep(list(1:m), 2)

                                        # here's the phylogeny (forget
                                        # the species names)
plot(phy)

                                        # and here's the associated
                                        # covariance matrix
image(as(Vphy, "sparseMatrix"))

                                        # put the covariance matrix in
                                        # a list, for model-input
                                        # purposes -- the idea is that
                                        # there might be other
                                        # covariance matrix (e.g. a
                                        # spatial one say).  it is
                                        # important that the list
                                        # element gets the name
                                        # `species` because this is
                                        # the name of the grouping
                                        # factor used in the model
                                        # formula below
covList <- list(species = Vphy)

                                        # formula interface! this
                                        # model has a fixed
                                        # interaction between the
                                        # environment and the trait
                                        # (with intercept and main
                                        # effects too), a random
                                        # environmental slope and
                                        # intercept with
                                        # (phylogenetic) correlations
                                        # across species.
form <- y ~ x * z + (x | species)

                                        # combine this formula with
                                        # the data and covariance
                                        # matrix to construct the
                                        # matrices and structures
                                        # required to create a
                                        # deviance function for
                                        # estimating the parameters of
                                        # the associated mixed model
parsedForm <- levelsCovFormula(form, df, covList = covList)

                                        # set the covariance
                                        # parameters to something more
                                        # interesting (i.e. with a
                                        # covariance between slope and
                                        # intercept)
covarSim <- c(0.5, -0.2, 0.5)
parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(covarSim))

                                        # update simulations to
                                        # reflect the new structure
X <- model.matrix(nobars(form), df) # fixed effects design matrix
Z <- t(parsedForm$Lambdat %*% parsedForm$Zt) # random effects design
                                             # matrix with
                                             # phylogenetic
                                             # covariances
fixefSim <- rnorm(ncol(X)) # fixed effects
u <- rnorm(ncol(Z)) # whitened random effects
p <- plogis(as.numeric(X %*% fixefSim + Z %*% u)) # probability of observation
dl$y <- rbinom(nrow(df), 1, p) # presence-absence data
df <- as.data.frame(dl) # reconstruct the data frame with new
                        # structured response

                                        # now look at the interesting
                                        # structure.  here's the
                                        # cholesky factor of the
                                        # species covariance
image(parsedForm$Lambdat)

                                        # and the species covariance
                                        # matrix itself.  The big four
                                        # blocks represent the 2-by-2
                                        # covariance between intercept
                                        # and slope.  The covariances
                                        # within these blocks
                                        # represent phylogenetic
                                        # covariance.  the pattern
                                        # here is more closely related
                                        # species have more similar
                                        # intercepts and slopes (red
                                        # blocks on the diagonal) but
                                        # more closely related species
                                        # also have stronger negative
                                        # correlations between slope
                                        # and intercept (blue blocks
                                        # on off diagonal)
image(crossprod(parsedForm$Lambdat))

                                        # here's the transposed random
                                        # effects model matrix.  those
                                        # are 1's for the intercepts
                                        # in the first 30 rows and the
                                        # environmental variable in
                                        # the second 30
image(parsedForm$Zt)

                                        # here's the full covariance
                                        # matrix (the large scale
                                        # blocks reflect phylogenetic
                                        # correlations and the
                                        # patterns within each block
                                        # are due to the environmental
                                        # variable)
image(fullCov <- t(parsedForm$Zt) %*% crossprod(parsedForm$Lambdat) %*% parsedForm$Zt)

                                        # here's a closeup of one of
                                        # the blocks
image(fullCov[1:10, 1:10])

                                        # a potential problem is that
                                        # this block is singular!
eigen(fullCov[1:10, 1:10])$values

                                        # in fact the rank of the full
                                        # 300 by 300 matrix is only 60
                                        # = 30 species times 2 column
                                        # model matrix
rankMatrix(fullCov)[1]

                                        # but then again so is the
                                        # standard non-phylogenetic
                                        # glmer model
gm <- glmer(form, df, binomial)
with(getME(gm, c("Zt", "Lambdat")), {
    covMatGm <- t(Zt) %*% crossprod(Lambdat) %*% Zt
    print(rankMatrix(covMatGm)[1])
    dim(covMatGm)
})

                                        # distribution of underlying
                                        # probabilities of occurrence
hist(p)

                                        # observed occurrence pattern
color2D.matplot(dl$y, xlab = "species", ylab = "sites", main = "abundance")

                                        # initial parameters
parInds <- list(covar = 1:3, fixef = 4:7, loads = NULL)
initPars <- c(covar = c(1, 0, 1), fixef = rep(0, 4))

                                        # make the deviance function
dfun <- mkGeneralGlmerDevfun(df$y, parsedForm$X,
                             parsedForm$Zt, parsedForm$Lambdat,
                             rep(1, nrow(df)), rep(0, nrow(df)),
                             initPars, parInds,
                             parsedForm$mapToCovFact, function(loads) NULL)

                                        # try it out
dfun(initPars)

                                        # optimize
opt <- bobyqa(initPars, dfun, lower = c(0, -Inf, 0, rep(-Inf, 4)),
              control = list(iprint = 4L))
names(opt$par) <- names(initPars)

                                        # compare with truth
cbind(estimated = opt$par, # estimated parameters
      true = c(covar = covarSim, fixef = fixefSim)) # true parameters
