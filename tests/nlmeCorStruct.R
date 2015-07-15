library(lme4ord)
library(ape)
library(multitable)

set.seed(1)

                                        # basic simulation design:
nSites <- 800
nSpec <- 10
form <- y ~ x + nlmeCorStruct(1 | species, corObj = corObj, sig = 1)
dataList <- dims_to_vars(data.list(y = 1 * (matrix(rnorm(nSites * nSpec), nSites, nSpec) > 0),
                                   x = rnorm(nSites), # environmental variable, x
                                   dimids = c("sites", "species")))
phy <- compute.brlen(rtree(n = nSpec), method = "Grafen", power = 0.8)
phy$tip.label <- dimnames(dataList)$species
corObj <- Initialize(corBrownian(1, phy), dropdl(dataList[1, ]))
image(repSparseCorFactor(corObj, sig = 1))
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

(gm <- strucGlmer(form, dataFrame, binomial, list(corObj = corObj)))
beta

