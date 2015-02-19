##' Parse mixed model formula for levels-covariance model
##'
##' @param formula mixed model formula
##' @param data data
##' @param family family
##' @param covList list of covariance matrices for random effects
##' grouping factors in \code{formula}
##' @param ... not used yet
##' @export
glmercFormula <- function(formula, data = NULL, family = binomial, covList, ...) {

                                        # get model matrix, grouping
                                        # factor, and term name
    reTrmsList <- lapply(findbars(formula), getModMatAndGrpFac, fr = data)
    names(reTrmsList) <- sapply(reTrmsList, "[[", "grpName")

                                        # append the covariance
                                        # matrices over the levels
                                        # associated with each random
                                        # effects term
    for(i in seq_along(reTrmsList)) {
        reTrmsList[[i]]$covMat <- covList[[reTrmsList[[i]]$grpName]]
        if(is.null(reTrmsList[[i]]$covMat)) { # if no cov mat, use
                                              # identity
            nli <- nlevels(reTrmsList[[i]]$grpFac)
            reTrmsList[[i]]$covMat <- diag(1, nli, nli)
        }
    }

                                        # extract the model matrices,
                                        # grouping factors, and cov
                                        # mats
    modMats <- lapply(reTrmsList, "[[", "modMat")
    grpFacs <- lapply(reTrmsList, "[[", "grpFac")
    covMats <- lapply(reTrmsList, "[[", "covMat")

                                        # set up constant grouping
                                        # factors and 1 by 1 cov mats
                                        # for mkTemplateReTrms (this
                                        # is where we could add
                                        # nesting)
    grpFacConst <- rep(list(as.factor(rep("const", nrow(data)))), length(grpFacs))
    covMatConst <- rep(list(matrix(1, 1, 1)), length(grpFacs))

                                        # construct the full random
                                        # effects model matrix,
                                        # relative covariance factor,
                                        # etc.
    re <- mkTemplateReTrms(modMats,
                           grpFacs, grpFacConst,
                           covMats, covMatConst)

                                        # organize output value
    return(c(re,
             list(X = model.matrix(nobars(formula), data),
                  y = model.response(model.frame(nobars(formula), data)),
                  flist = simplifyFacList(grpFacs))))
}
