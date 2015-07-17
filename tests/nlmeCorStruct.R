library(lme4ord)
library(ape)
library(nlme)

set.seed(1)

                                        # simulation design
nLakes   <- 800
nSpecies <- 10
lakes   <- data.frame(lakes   = paste("lake",    1:nLakes,   sep = ""))
species <- data.frame(species = paste("species", 1:nSpecies, sep = ""))
rownames(species) <- species$species
dataFrame <- merge(species, lakes)
dataFrame$envVar <- rnorm(nLakes)[as.numeric(dataFrame$lakes)]
phy <- compute.brlen(rtree(n = nSpecies), method = "Grafen", power = 0.8)
phy$tip.label <- as.character(unique(dataFrame$species))
corObj <- Initialize(corBrownian(1, phy), species)


                                        # model design
form <- respVar ~ envVar + nlmeCorStruct(1 | species, corObj = corObj, sig = 1)
pform <- strucParseFormula(form, dataFrame, addArgs = list(corObj = corObj))
pform$initPars[2:3] <- beta <- c(-1, 2)

                                        # simulate response
respVar <- simulate(pform, nsim = 1, seed = 1,
                    weights = rep(1, nrow(dataFrame)),
                    family = binomial())
dataFrame$respVar <- as.numeric(respVar)

                                        # fit model
(gm <- strucGlmer(form, dataFrame, binomial, list(corObj = corObj)))

                                        # plot covariance and phylogeny
image(crossprod(getReTrm(gm, "species.nlmeCorStruct")$Lambdat))
plot(phy)

