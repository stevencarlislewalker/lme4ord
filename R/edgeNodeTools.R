##' Relationships between tips and edges
##'
##' \code{findFromNode} is a recurrsive function for computing the
##' path from a node back in time through its ancestor nodes.
##' \code{findEdgesFromPath} returns a \code{logical} vector
##' indicating the edges associated with a particular path.
##' \code{edgeTipIndicator} returns a \code{logical} matrix giving the
##' relationship between tips and the edges associated with their
##' history.
##' 
##' @param to vector of node indices
##' @param edge phylogenetic edge matrix
##' @param path vector giving the path from a node back in time
##' through its ancestor nodes (output of \code{findFromNode})
##' @param object either a \code{phylo} or \code{matrix} object
##' @rdname edgeTipIndicator
##' @export
##' @examples
##' set.seed(1)
##' phy <- rtree(5)
##' plot(phy)
##' edgelabels()
##' edgeTipIndicator(phy)
findFromNode <- function(to, edge) {
    lastTo <- to[length(to)]
    newTo <- edge[edge[, 2] == lastTo, ][1]
    if(is.na(newTo)) return(to)
    findFromNode(c(to, newTo), edge)
}
##' @rdname edgeTipIndicator
##' @export
findEdgesFromPath <- function(path, edge) {
    col1 <- edge[, 1] %in% path[-1           ]
    col2 <- edge[, 2] %in% path[-length(path)]
    col1 & col2
}
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator <- function(object, ...) {
    UseMethod("edgeTipIndicator")
}
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.default <- function(object, ...) {
    edgeTipIndicator(as.matrix(object))
}
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.matrix <- function(object, ntip, ...) {
    if(ncol(object) != 2L) stop("not an edge matrix")
    sapply(lapply(1:ntip, findFromNode, object),
           findEdgesFromPath, object)
}
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.phylo <- function(object, ...) {
    ans <- edgeTipIndicator(object$edge, Ntip(object))
    colnames(ans) <- object$tip.label
    return(ans)
}
