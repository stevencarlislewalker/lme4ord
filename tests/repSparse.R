library(lme4ord)

# ------------------------------------------------------------
# construction
# ------------------------------------------------------------

X <- strucMatrixTri(c(1.5, 1), -0.1)
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
# subset
# ------------------------------------------------------------

set.seed(1)
n <- 8; m <- 5
X1 <- strucMatrix(1:m,
               rep(1:2, c(2, m - 2)),
               1:m, rep(1, m))
fac <- factor(letters[rep(sample(m), n)])
levels(fac) <- levels(fac)[sample(m)]
Z <- subset(X1, as.numeric(fac))
image(update(Z, rnorm(m)))
image(update(Z, rnorm(m)))

# ------------------------------------------------------------
# updates
# ------------------------------------------------------------

as.matrix(update(X, c(1, 1, -0.2)), TRUE)
as.matrix(update(X, diagVals = c(1.2, 1.5), offDiagVals = 10), TRUE)
as.matrix(update(Xrep, diagVals = c(1.2, 1.5), offDiagVals = 10), TRUE)

update(XkronY, c(rnorm(5), rnorm(3)))


# ------------------------------------------------------------
# specials
# ------------------------------------------------------------

X <- strucMatrixVarWithCovariate(rnorm(5), rnorm(20), gl(5, 4),
                               mkVarPowerTrans)


                                        # compound symmetry
X <- strucMatrixCompSymm(1, -0.2, 5)
                                        # check that cholesky follows
                                        # as.matrix is equivalent to
                                        # as.matrix follows cholesky
stopifnot(all.equal(chol(as.matrix(X, TRUE)),
                    t(as.matrix(chol(X), TRUE))))
                                        # ditto for the update method
stopifnot(all.equal(chol(as.matrix(update(X, c(1.2, -0.2)), TRUE)),
                    t(as.matrix(update(chol(X), c(1.2, -0.2)), TRUE))))
                                        # and if we take the transpose
                                        # of the chol(strucMatrix)
                                        # itself
stopifnot(all.equal(chol(as.matrix(update(X, c(1.2, -0.2)), TRUE)),
                    as.matrix(update(t(chol(X)), c(1.2, -0.2)), TRUE)))
                                        # check if chol updating works
                                        # with kron
Y <- strucMatrixTri(c(1, 2), -0.1)
stopifnot(all.equal(update(kron(Y, chol(X)), c(2, -0.1, 4, 1, 3)),
                    kron(update(Y, c(4, 1, 3)), update(chol(X), c(2, -0.1)))))




