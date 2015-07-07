library(lme4ord)

glmer(incidence/size ~ period + (1 | herd), cbpp, binomial, weights = cbpp$size)
strucGlmer(prop ~ period + lme4(1 | herd),
           within(cbpp, prop <- incidence/size),
           binomial, weights = cbpp$size)

strucGlmer(prop ~ period + lme4(1 | herd),
           within(cbpp, prop <- incidence/size),
           binomial, weights = cbpp$size,
           penCovar = function(x) 5 * abs(x))


