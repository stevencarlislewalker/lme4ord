library(lme4ord)
#library(Matrix)
#library(vegan)

set.seed(1)

                                        # simulation design
nLakes   <- 50
nSpecies <- 20
lakes   <- paste("lake",    1:nLakes,   sep = "")
species <- paste("species", 1:nSpecies, sep = "")
dataFrame <- expand.grid(lakes   = lakes,
                         species = species)

                                        # model design
form  <- respVar ~ 1 + lme4(1 | lakes) + lme4(1 | species)
form2 <- respVar ~ 1 + lme4(1 | lakes)

                                        # simulations
pform <- strucParseFormula(form, dataFrame)
dataFrame$respVar <- simulate(pform, nsim = 1, seed = 1,
                              weights = rep(1, nrow(dataFrame)),
                              family = binomial())

gm <- strucGlmer(form, dataFrame, binomial)
gmu <- update(gm, . ~ . -lme4(1 | species), data = dataFrame)
gm2 <- strucGlmer(form2, dataFrame, binomial)
stopifnot(all(pars(gmu) == pars(gm2)))
