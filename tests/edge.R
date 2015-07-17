library(lme4ord)
library(multitable)
library(ape)

set.seed(1)

                                        # basic simulation design:
nSites <- 10
nSpec <- 100
form <- respVar ~ envVar * trait +        # environment by trait interaction
    (0 + trait | sites) +                 # trait effects vary over sites
    edge(1     | species, phy = phy) +    # the intercept is phylogenetically correlated
    edge(0 + envVar | species, phy = phy) # the environment effect is phylogenetically correlated
dataList <- dims_to_vars(data.list(respVar = 1*(matrix(rnorm(nSites*nSpec), nSites, nSpec) > 0),
                                   envVar = rnorm(nSites), # environmental variable
                                   trait  = rnorm(nSpec),  # trait
                                   dimids = c("sites", "species")))
phy <- compute.brlen(rtree(n = nSpec), method = "Grafen", power = 0.8)
phy$tip.label <- dimnames(dataList)$species
dataFrame <- as.data.frame(dataList)

pform <- strucParseFormula(form, dataFrame, addArgs = list(phy = phy))
pform$initPars[pform$parInds$fixef] <- beta <- rnorm(ncol(pform$fixed))

respVar <- simulate(pform, nsim = 1, seed = 1,
                    weights = rep(1, nrow(dataFrame)),
                    family = binomial())
dataFrame$respVar <- c(respVar)

                                        # fit model
(gm <- strucGlmer(form, dataFrame, binomial, list(phy = phy)))
beta

