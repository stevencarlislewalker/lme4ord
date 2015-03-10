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
##' (ee <- edgeTipIndicator(phy))
##' if (require(Matrix)) {
##'    image(Matrix(ee),sub="",xlab="tips",ylab="branches")
##' }
findFromNode <- function(to, edge) {
    lastTo <- to[length(to)]
    newTo <- edge[edge[, 2] == lastTo, ][1]
    if(is.na(newTo)) return(to)
    findFromNode(c(to, newTo), edge)
}
##' @rdname edgeTipIndicator
##' @param scale scaling factor for edges (i.e., vector of branch lengths)
##' @export
findEdgesFromPath <- function(path, edge, scale=1) {
    col1 <- edge[, 1] %in% path[-1           ]
    col2 <- edge[, 2] %in% path[-length(path)]
    (col1 & col2)*scale
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
edgeTipIndicator.matrix <- function(object, ntip, scale=1, ...) {
    if(ncol(object) != 2L) stop("not an edge matrix")
    sapply(lapply(1:ntip, findFromNode, object),
           findEdgesFromPath, object, scale)
}
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.phylo <- function(object, ...) {
    ## FIXME: do all phylo objects have edge lengths or can they
    ##  be implicit?
    ans <- edgeTipIndicator(object$edge, Ntip(object), object$edge.length)
    colnames(ans) <- object$tip.label
    return(ans)
}
