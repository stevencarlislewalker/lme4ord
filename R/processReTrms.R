##' @param modMat raw random effects model matrix
##' @param grpFac grouping factor
##' @param grpName grouping factor name
##' @rdname reStruct
##' @export
reStructIdentity <- function() {
    list(getReTrms =
         function(modMat, grpFac, grpName) {
             Zt <- kr(t(as.repSparse(modMat)), as.repSparse(grpFac))
             return(list(Zt = resetTransConst(Zt),
                         Lambdat = repSparseIdent(nlevels(grpFac))))
         },
         printReTrms =
         function() {
             print("identity reStruct")
         })
}
