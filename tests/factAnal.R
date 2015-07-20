library(lme4ord)
library(Matrix)
library(vegan)

set.seed(1)

                                        # simulation design
nLakes   <- 50
nSpecies <- 30
lakes   <- paste("lake",    1:nLakes,   sep = "")
species <- paste("species", 1:nSpecies, sep = "")
dataFrame <- expand.grid(lakes   = lakes,
                         species = species)

                                        # model design
form <- respVar ~ 1 + 
    (1 | lakes) +
    (1 | species) +
    factAnal(0 + lakes | species, nAxes = 2, seed = 1)

                                        # simulations 
pform <- strucParseFormula(form, dataFrame)
respVar <- simulate(pform, nsim = 1, seed = 1,
                    weights = rep(1, nrow(dataFrame)),
                    family = binomial())
dataFrame$respVar <- as.numeric(respVar)
respMat <- matrix(dataFrame$respVar,
                  nLakes, nSpecies)
print(dataFrame)
image(as(respMat, "sparseMatrix"))

                                        # regular GLMM
gm0 <- glmer(noSpecials(form), dataFrame, binomial)
print(gm0)

                                        # structured GLMM with
                                        # factor analysis
gm <- strucGlmer(form, dataFrame, binomial, optMaxit = 50000,
                 penLoads = mkPenLpNorm(p = 2, lambda = 1))
print(gm)

fixef(gm) # 'true' value 0
covarPerTerm(gm) # 'true' values both 1

                                        # extract factor analysis
gmFA <- getReTrm(gm, "species.factAnal")

                                        # extract scores and
                                        # loadings
gmScores <- scores(gmFA)
head(loadings(gmScores))
head(factors(gmScores))

                                        # look at biplot
biplot(gmFA)

                                        # test the fit of model
fitMat <- matrix(plogis(fitted(gm)), nLakes, nSpecies)
j <- 14
boxplot(fitMat[,j] ~ respMat[,j],
        horizontal = TRUE, ylim = c(0, 1))

                                        # reordering the matrix
ordMat <- matrix(respMat,
                 nLakes,
                 nSpecies)[order(factors(gmScores)[,2]),
                           order(loadings(gmScores)[,2])]
image(as(respMat, "sparseMatrix"))
image(as(ordMat, "sparseMatrix"))
