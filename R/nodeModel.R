##' Tip-node bigraph adjacency matrix
##'
##' @param phy phylo object
##' @return adjacency matrix
##' @importFrom reshape2 acast melt
##' @importFrom ape Ntip Nnode prop.part
##' @export
tipNodeRelationship <- function(phy) {
    nodeNms <- if(is.null(nl <- phy$node.label)) Ntip(phy) + (1:Nnode(phy)) else nl
    propPartPhy <- setNames(unclass(prop.part(phy)), nodeNms)
    return(1 * (acast(melt(propPartPhy), value ~ L1, fill = 0) > 0L))
}
