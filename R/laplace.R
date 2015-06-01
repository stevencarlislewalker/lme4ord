##' Laplace approximation
##'
##' Compute the Laplace approximation of the log-likelihood of a
##' generalized linear mixed model.
##'
##' @param object fitted model object
##' @param ... not used
##' @export
laplace <- function(object, ...) {
    UseMethod("laplace")
}

##' @rdname laplace
##' @export
laplace.default <- function(object, ...) {
    resp <- object$resp
    pp <- object$pp
    if(is.null(resp) || is.null(pp)) stop("object must contain resp or pp")
    if(!inherits(resp, "glmResp") || !inherits(pp, "merPredD")) {
        stop("resp must be a glmResp object and ",
             "pp must be a merPredD object")
    }
    -0.5 * (resp$aic() + pp$sqrL(1) + pp$ldL2())
}

##' @rdname laplace
##' @export
laplace.glmerMod <- function(object, ...) {
    laplace(list(pp = object@pp, resp = object@resp))
}

##' @rdname laplace
##' @export
laplace.strucGlmer <- function(object, ...) {
    laplace(object$parsedForm$devfunEnv)
}
