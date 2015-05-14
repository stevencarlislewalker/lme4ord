triKron <- function(A, B) {
    n <- dim(A)[1]
    m <- dim(B)[1]
    sparseMatrix(i = 1 + as.integer(outer(B@i, m * A@i, "+")),
                 j = 1 + as.integer(outer(B@j, m * A@j, "+")),
                 x = as.numeric(outer(B@x, A@x, "*")),
                 dims = rep(n * m, 2),
                 giveCsparse = FALSE)
}

triDsum <- function(A, B) {
    n <- dim(A)[1]
    m <- dim(B)[1]
    sparseMatrix(i = 1 + c(A@i, n + B@i),
                 j = 1 + c(A@j, n + B@j),
                 x = c(A@x, B@x),
                 dims = rep(n + m, 2),
                 giveCsparse = FALSE)
}

library(Matrix)

R1 <- as(t(chol(crossprod(rmat(4, 3)))), "TsparseMatrix")
T1 <- as(t(chol(crossprod(rmat(3, 2)))), "TsparseMatrix")
S1 <- as(t(chol(crossprod(rmat(5, 4)))), "TsparseMatrix")
R2 <- as(t(chol(crossprod(rmat(4, 3)))), "TsparseMatrix")
T2 <- as(t(chol(crossprod(rmat(3, 2)))), "TsparseMatrix")
S2 <- as(t(chol(crossprod(rmat(5, 4)))), "TsparseMatrix")

image(triDsum(Reduce(triKron, list(R1, T1, S1)),
              Reduce(triKron, list(R2, T2, S2))))




image(R %x% T %x% S)
image(Reduce(triKron, list(R, T, S)))

RxT <- sparseMatrix(i = 1 + as.integer(outer(T@i, dim(T)[1] * R@i, "+")),
                    j = 1 + as.integer(outer(T@j, dim(T)[1] * R@j, "+")),
                    x = as.numeric(outer(T@x, R@x, "*")),
                    dims = rep(dim(R)[1] * dim(T)[1], 2),
                    giveCsparse = FALSE)


image(triKron(triKron(R, T), S))
image(R %x% T %x% S)

image(triDsum(triKron(triKron(R, T), S), triKron(triKron(R, T), S)))



cbind(RxT@i, (R%x%T)@i)
cbind(RxT@j, (R%x%T)@j)
cbind(RxT@x, (R%x%T)@x)
image(t(RxT))



