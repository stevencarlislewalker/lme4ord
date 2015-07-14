library(lme4ord)
library(multitable)
library(vegan)

data(mite)
data(mite.env)
data(mite.xy)

dl <- dims_to_vars(data.list(mite = as.matrix(mite), mite.env, mite.xy,
                             dimids = c("sites", "species")))
dl <- aperm(dl, c(2, 1)) ## better sparsity properties if species
                         ## dimension comes first
dl$sites <- factor(dl$sites, dl$sites)
summary(dl)
df <- as.data.frame(dl)

miteDist <- dist(cbind(dl$x, dl$y))

form <- mite ~ 1 + (1 | species) + 
    expDecay(1 | sites, distMat = distMat,
             minCov = 1e-3, distCutoff = 2)

(gm <- strucGlmer(form, df, poisson, addArgs = list(distMat = miteDist)))

par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
with(covExpDecay(covarPerTerm(gm)$sites.expDecay, distCutoff = 2),
     plot(edgeDists, edgeCovs, type = "l", xlim = c(0, max(miteDist))))
segments(2, 0, max(miteDist), 0)
hist(miteDist, 30)

