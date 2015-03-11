##' Make edge-model random effects term
##'
##' @param modMat a model matrix
##' @param grpFac grouping factor (typically dimension IDs of a melted
##' response matrix)
##' @param edgeMat matrix giving the relationships between the levels
##' of \code{grpFac} and the edges of the tree (the levels should be
##' the tips of the tree)
##' @param giveCsparse use compressed-form \code{CsparseMatrix}
##' representations (\code{TRUE}), or triplet-form
##' \code{TsparseMatrix} representations (\code{FALSE})?
##' @export
mkTemplateReTrm <- function(modMat, grpFac, edgeMat,
                            giveCsparse = TRUE) {

    matClass <- if(giveCsparse) "CsparseMatrix" else "TsparseMatrix"

    J <- as(as(grpFac, "sparseMatrix"), matClass)
    Zt <- KhatriRao(edgeMat %*% J, t(modMat))
    

}
