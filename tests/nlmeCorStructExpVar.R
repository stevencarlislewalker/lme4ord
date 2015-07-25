library(lme4ord)
library(lme4)

if(require(nlme)) {

set.seed(1)


                                        # simulation design
nLakes   <- 20
nSpecies <- 10
lakes   <- sprintf("lake%02d",    1:nLakes)
species <- sprintf("species%02d", 1:nSpecies)
dataFrame <- expand.grid(lakes   = lakes,
                         species = species)
spatDat <- data.frame(lakes = lakes, x = rnorm(nLakes), y = rnorm(nLakes))


corObj <- Initialize(corExp(1, form = ~ x + y), data = spatDat)
coef(corObj, FALSE) #   constrained form
coef(corObj, TRUE)  # unconstrained form

coef(corObj, TRUE)[] <- 0.5

form <- respVar ~ 1 +
    (1 | species) +
    nlmeCorStruct(1 | lakes, corObj = corObj)

pform <- strucParseFormula(form, dataFrame, addArgs = list(corObj = corObj))

respVar <- simulate(pform, nsim = 1, seed = 1,
                    weights = rep(1, nrow(dataFrame)),
                    family = binomial())

dataFrame$respVar <- as.numeric(respVar)

gm <- strucGlmer(form, dataFrame, binomial,
                 addArgs = list(corObj = corObj),
                 optMaxit = 50000)

image(as(VarCorr(getReTrm(gm, "lakes.nlmeCorStruct")), "sparseMatrix"))

specDat <- data.frame(species = species, specRanef = ranef(gm)$species)
spatDat$lakeRanef <- ranef(gm)$lakes

if(require(ggplot2)) {
    ggplot(merge(merge(dataFrame, spatDat), specDat), aes(x, y)) +
        geom_point(aes(colour = as.factor(respVar), size = lakeRanef), alpha = 0.8) +
            facet_wrap( ~ species, ncol = 2) +
                theme_bw()
    
    ggplot(merge(merge(dataFrame, spatDat), specDat)) +
        geom_point(aes(x, y, size = lakeRanef))
}

}
