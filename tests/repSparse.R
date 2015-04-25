library(lme4ord)

# ------------------------------------------------------------
# construction
# ------------------------------------------------------------

X <- repSparseTri(c(1.5, 1), -0.1)
Xrep <- rep(X, 3, type = "diag")
as.matrix(X, TRUE)
as.matrix(Xrep, TRUE)

set.seed(1)
Y <- rRepSparse(2, 7, 5, 9)


# ------------------------------------------------------------
# binds and reps
# ------------------------------------------------------------


# ------------------------------------------------------------
# products
# ------------------------------------------------------------

XkronY <- kron(X, Y)

# ------------------------------------------------------------
# updates
# ------------------------------------------------------------

as.matrix(update(X, c(1, 1, -0.2)), TRUE)
as.matrix(update(X, diagVals = c(1.2, 1.5), offDiagVals = 10), TRUE)
as.matrix(update(Xrep, diagVals = c(1.2, 1.5), offDiagVals = 10), TRUE)

XkronY
update(XkronY, c(rnorm(5), rnorm(3)))


# ------------------------------------------------------------
# specials
# ------------------------------------------------------------

X <- repSparseVarWithCovariate(rnorm(5), rnorm(20), gl(5, 4),
                               mkVarPowerTrans)
image(X)

pr <- rnorm(5)
upX <- update(X, pr)
image(upX)
with(environment(upX$trans), {
    cbind(covariate,
          pr[as.numeric(grpFac)])
})


                                        # compound symmetry
X <- repSparseCompSymm(1, -0.2, 5)
                                        # check that cholesky follows
                                        # as.matrix is equivalent to
                                        # as.matrix follows cholesky
stopifnot(all.equal(as(chol(as.matrix(X, TRUE)), "dgCMatrix"),
                    as.matrix(t(chol(X)), TRUE)))
                                        # ditto for the update method
stopifnot(all.equal(as(chol(as.matrix(update(X, c(1.2, -0.2)), TRUE)), "dgCMatrix"),
                    t(as.matrix(update(chol(X), c(1.2, -0.2)), TRUE))))
                                        # and if we take the transpose
                                        # of the chol(repSparse)
                                        # itself
stopifnot(all.equal(as(chol(as.matrix(update(X, c(1.2, -0.2)), TRUE)), "dgCMatrix"),
                    as.matrix(update(t(chol(X)), c(1.2, -0.2)), TRUE)))





