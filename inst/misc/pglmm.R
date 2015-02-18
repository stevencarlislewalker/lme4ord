library("lme4ord")
library("minqa")
library("plotrix")
library("ape")

                                        # initial simulation (will
                                        # change y below to have
                                        # structure)
set.seed(10)
n <- 10
m <- 30
dl <- dims_to_vars(data.list(y = 1 * (rmat(n, m) > 0),
                             x = rnorm(n), z = rnorm(m),
                             dimids = c("sites", "species")))
df <- as.data.frame(dl)
head(df)

                                        # make up some silly
                                        # covariance matrix over the
                                        # species
phy <- rtree(n = m)
phy <- compute.brlen(phy, method = "Grafen", power = 0.5)
Vphy <- stanCov(vcv(phy))
covList <- list(species = Vphy)
plot(phy)

                                        # formula interface!
parsedForm <- levelsCovFormula(y ~ x + (x | species), df, covList = covList)

                                        # set the covariance
                                        # parameters to something more
                                        # interesting (i.e. with a
                                        # covariance between slope and
                                        # intercept)
thetaSim <- c(0.5, -0.2, 0.5)
parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(thetaSim))

                                        # here's the cholesky factor
                                        # of the species covariance
image(parsedForm$Lambdat)

                                        # here's the transposed random
                                        # effects model matrix
image(parsedForm$Zt)

                                        # update simulations to
                                        # reflect the new structure
X <- model.matrix(y ~ x, df)
Z <- t(parsedForm$Lambdat %*% parsedForm$Zt)
beta <- rnorm(ncol(X))
u <- rnorm(ncol(Z))
p <- plogis(as.numeric(X %*% beta + Z %*% u))
dl$y <- rbinom(nrow(df), 1, p)
df <- as.data.frame(dl)

                                        # distribution of underlying
                                        # probabilities of occurrence
hist(p)

                                        # observed occurrence pattern
color2D.matplot(dl$y, xlab = "species", ylab = "sites", main = "abundance")

                                        # initial parameters
parInds <- list(covar = 1:3, fixef = 4:5, loads = NULL)
initPars <- c(covar = c(1, 0, 1), fixef = c(0, 0))

                                        # make the deviance function
dfun <- mkGeneralGlmerDevfun(df$y, parsedForm$X,
                             parsedForm$Zt, parsedForm$Lambdat,
                             rep(1, nrow(df)), rep(0, nrow(df)),
                             initPars, parInds,
                             parsedForm$mapToCovFact, function(loads) NULL)

                                        # try it out
dfun(initPars)

                                        # optimize
opt <- bobyqa(initPars, dfun, lower = c(0, -Inf, 0, -Inf, -Inf),
              control = list(iprint = 4L))
names(opt$par) <- names(initPars)

                                        # compare with truth -- works
                                        # ok, but certainly not all
                                        # the time -- e.g. set.seed(1)
opt$par
c(covar = thetaSim, fixef = beta)
