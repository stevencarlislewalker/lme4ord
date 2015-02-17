##' Factors with covariance matrices over the levels
##'
##' @param fac a factor
##' @param covMat a covariance matrix over the levels (if missing use
##' identity matrix)
##' @param stan standardize covariance matrix to determinant one?
##' @return a \code{factorCov} object
##' @rdname factorCov
##' @export
factorCov <- function(fac, covMat, stan = TRUE) {
    fac <- as.factor(fac)
    if(missing(covMat)) {
        covMat <- diag(1, nlevels(fac), nlevels(fac))
    } else if(stan) {
        covMat <- stanCov(covMat)
    }
    if(!(nlevels(fac) == nrow(covMat))) stop("covMat must be over levels(fac)")
    attr(fac, "covMat") <- covMat
    class(fac) <- c("factorCov", "factor")
    return(fac)
}
##' @rdname factorCov
##' @export
print.factorCov <- function(x, ...) {
    ret <- x
    attr(x, "covMat") <- NULL
    class(x) <- "factor"
    print(as.factor(x), ...)
    cat("Note: this factor includes a covariance matrix over its levels, attr(., \"covMat\")\n")
}
##' Standardize covariance matrix to determinant one
##'
##' @param covMat
##' @export
stanCov <- function(covMat) {
    covMat / (det(covMat)^(1/nrow(covMat)))
}
