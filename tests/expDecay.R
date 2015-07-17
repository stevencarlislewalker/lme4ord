library(lme4ord)
library(lme4)

                                        # simulation design
nLakes   <- 50
nSpecies <- 30
lakes   <- paste("lake",    1:nLakes,   sep = "")
species <- paste("species", 1:nSpecies, sep = "")
dataFrame <- expand.grid(lakes   = lakes,
                         species = species)

                                        # model design
form <- respVar ~ 1 + 
    (1 | species) +
    expDecay(1 | lakes, distMat = distMat, minCov = 1e-3, distCutoff = 1.5)
distMat <- dist(matrix(rnorm(nLakes * 2), nLakes, 2))

                                        # simulations
pform <- strucParseFormula(form, dataFrame, addArgs = list(distMat = distMat))
respVar <- simulate(pform, nsim = 1, seed = 1,
                    weights = rep(1, nrow(dataFrame)),
                    family = binomial())
dataFrame$respVar <- as.numeric(respVar)
respMat <- matrix(dataFrame$respVar,
                  nLakes, nSpecies)

                                        # structured GLMM with
                                        # factor analysis
gm <- strucGlmer(form, dataFrame, binomial,
                 addArgs = list(distMat = distMat),
                 optMaxit = 50000)
print(gm)

par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
with(covExpDecay(covarPerTerm(gm)$lakes.expDecay, distCutoff = 1.5),
     plot(edgeDists, edgeCovs, type = "l", xlim = c(0, max(distMat))))
segments(1.5, 0, max(distMat), 0)
hist(distMat, 30)

