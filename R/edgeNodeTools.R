##' Relationships between tips and edges
##'
##' \code{findFromNode} is a recurrsive function for computing the
##' path from a node back in time through its ancestor nodes.
##' \code{findEdgesFromPath} returns a \code{logical} vector
##' indicating the edges associated with a particular path.
##' \code{edgeTipIndicator} returns a \code{logical} matrix giving the
##' relationship between tips and the edges associated with their
##' history.  \code{reorderPhylo} is a convenience function for
##' ordering the edges of a \code{phylo} object in a way that makes it
##' easier to build edge-based phylogenetic models.
##' 
##' @param node vector of node indices
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
##' (ee <- edgeTipIndicator(phy))
##' if (require(Matrix)) {
##'    image(Matrix(ee),sub="",xlab="tips",ylab="branches")
##' }
findPathFromNode <- function(node, edge) {
    lastNode <- node[length(node)]
    newNode <- edge[edge[, 2] == lastNode, ][1]
    if(is.na(newNode)) return(node)
    findPathFromNode(c(node, newNode), edge)
}

##' @rdname edgeTipIndicator
##' @param scale scaling factor for edges (i.e., vector of branch lengths)
##' @export
findEdgesFromPath <- function(path, edge, scale = 1) {
    col1 <- edge[, 1] %in% path[-1           ]
    col2 <- edge[, 2] %in% path[-length(path)]
    (col1 & col2)*scale
}

##' @param ... not used
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

##' @param ntip number of tips
##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.matrix <- function(object, ntip, scale = NULL, ...) {
    if(ncol(object) != 2L) stop("not an edge matrix")
    ans <- sapply(lapply(1:ntip, findPathFromNode, object),
                  findEdgesFromPath, object)
    attr(ans, "scale") <- scale
    return(ans)
}

##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.phylo <- function(object, ...) {
    scl <- if(is.null(object$edge.length)) NULL else object$edge.length
    ans <- edgeTipIndicator(object$edge, Ntip(object), scl)
    colnames(ans) <- object$tip.label
    return(ans)
}

##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.hclust <- function(object, ...) {
    edgeTipIndicator(as.phylo(object))
}

##' @rdname edgeTipIndicator
##' @export
edgeTipIndicator.mst <- function(object, ...) {
    sparseMat <- as(unclass(object), "TsparseMatrix")
    lowerInds <- sparseMat@i > sparseMat@j
    ans <- sparseMatrix(i = rep(seq_along(which(lowerInds)), 2),
                        j = c(sparseMat@i[lowerInds] + 1, sparseMat@j[lowerInds] + 1),
                        x = rep(1, 2 * sum(lowerInds)))
    colnames(ans) <- colnames(object)
    return(as.matrix(ans))
}

##' @rdname edgeTipIndicator
##' @seealso \code{\link{reorder.phylo}}
##' @export
reorderPhylo <- function(object, ...) {
    structure(within(unclass(object), {
        for(i in 2:1) {
            ord <- order(edge[, i], decreasing = TRUE)
            edge <- edge[ord, ]
            if(!is.null(object$edge.length)) {
                object$edge.length <- object$edge.length[ord]
            }
        }
    }), class = class(object))
}

