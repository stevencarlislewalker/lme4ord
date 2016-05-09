library(lme4ord)
stopifnot(all.equal(familySimFun(gaussian),
                    familySimFun(gaussian())))
stopifnot(all.equal(familySimFun("gaussian"),
                    familySimFun(gaussian())))

set.seed(1)
stopifnot(all.equal(familySimFun(binomial)(rep(1, 10), 10, runif(10)),
                    c(0, 0, 0, 1, 0, 1, 1, 0, 1, 0)))

set.seed(1)
stopifnot(all.equal(familySimFun(poisson)(rep(1, 10), 10, 1:10),
                    c(0L, 1L, 3L, 7L, 3L, 9L, 11L, 9L, 10L, 14L)))
