library(lme4ord)
library(multitable)
library(ape)

set.seed(1)

                                        # basic simulation design:
nSites <- 10
nSpec <- 100
form <- y ~ x * z +                    # environment by trait interaction
    (0 + z | sites) +                  # trait effects vary over sites
    edge(1     | species, phy = phy) + # the intercept is phylogenetically correlated
    edge(0 + x | species, phy = phy)   # the environment effect is phylogenetically correlated
dataList <- dims_to_vars(data.list(y = 1 * (matrix(rnorm(nSites * nSpec), nSites, nSpec) > 0),
                                   x = rnorm(nSites), # environmental variable, x
                                   z = rnorm(nSpec),  # trait, z
                                   dimids = c("sites", "species")))
phy <- compute.brlen(rtree(n = nSpec), method = "Grafen", power = 0.8)
phy$tip.label <- dimnames(dataList)$species
dataFrame <- as.data.frame(dataList)

                                        # simualate response variable:
parsedForm <- strucParseFormula(form, dataFrame, addArgs = list(phy = phy))
beta <- rnorm(ncol(parsedForm$fixed))
u <- rnorm(nrow(parsedForm$Lambdat))
dataFrame$y <- with(parsedForm, {
    fe <- as.numeric(fixed %*% beta)
    re <- as.numeric(t(Lambdat * Zt) %*% u)
    mu <- plogis(fe + re)
    rbinom(nSites * nSpec, 1, mu)
})

                                        # fit model
(gm <- strucGlmer(form, dataFrame, binomial, list(phy = phy)))
beta

