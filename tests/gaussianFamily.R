if(FALSE) {

library(lme4ord)

mLmer <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
mGlmer <- glmer(Reaction ~ Days + (Days | Subject), sleepstudy, family = freeDispGaussian())
mStrucGlmer <- strucGlmer(Reaction ~ Days + lme4(Days | Subject), sleepstudy, family = gaussian)

pForm <- glFormula(Reaction ~ Days + (Days | Subject), sleepstudy, family = freeDispGaussian())
dfun0 <- do.call(mkDglmerDevfun, pForm)
opt0 <- lme4:::optimizeDglmer(dfun0)
dfun <- updateGlmerDevfun(dfun0, pForm$reTrms)

opt <- lme4:::optimizeDglmer(dfun, stage = 2)



dfun1 <- do.call(mkGlmerDevfun, mGlmer)

optimizeGlmer(dfun1)

parsedForm <- strucParseFormula(Reaction ~ Days + lme4(Days | Subject), sleepstudy)

devf <- pirls(X = model.matrix(parsedForm),
              y = response(parsedForm),
              Zt = ranefModMat(parsedForm),
              Lambdat = relCovFact(parsedForm),
              thfun = parsedForm$mapToCovFact,
              theta = covar(parsedForm),
              weights = rep(1, nobs(parsedForm)),
              offset = rep(0, nobs(parsedForm)),
              eta = response(parsedForm),
              family = freeDispGaussian)

devf(c(pars(parsedForm), 1))

opt <- optim(c(pars(parsedForm), 1), devf, method = "L-BFGS-B")

mLmer
mGlmer
mStrucGlmer

getME(mLmer, "theta")
getME(mGlmer, "theta")
covar(mStrucGlmer)
opt$par[1:3]

mGlmer@resp$variance()

parLength(mStrucGlmer)["fixef"]

VarCorr(mStrucGlmer)$Subject.lme4

sqrt(sum(mStrucGlmer$parsedForm$devfunEnv$resp$devResid())/180)

sqrt(sum(residuals(mStrucGlmer)^2)/(nobs(mStrucGlmer)))

strucDims(mStrucGlmer)

}
