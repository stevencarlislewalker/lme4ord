##' Parse a mixed model formula
##'
##' @param formula mixed model formula
##' @param data an object coercible to data frame
##' @param ... additional parameters to \code{\link{as.data.frame}}
##' @export
generalParseFormula <- function(formula, data, ...) {
    data <- as.data.frame(data, ...)
                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

    return(list(response = model.response(model.frame(nobars(formula), data)),
                fixed    = model.matrix(nobars(formula), data),
                random   = reTrmsList))
}
