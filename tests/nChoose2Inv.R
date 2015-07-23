library(lme4ord)
xx <- floor(2^seq(0, 100, length = 100))
yy <- choose(xx, 2)
stopifnot(all.equal(xx, nChoose2Inv(yy)))

