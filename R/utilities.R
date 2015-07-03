##' Standardize covariance matrix to determinant one
##'
##' @param covMat covariance matrix
##' @export
stanCov <- function(covMat) {
    covMat / (det(covMat)^(1/nrow(covMat)))
}

##' Get model matrix and grouping factor
##'
##' This is kind of like \code{\link{model.matrix}} but for random
##' effects terms.
##' 
##' @param bar random effect language object (e.g. \code{x | g})
##' @param fr model frame
##' @return list with model matrix and grouping factor
##' @export
getModMatAndGrpFac <- function(bar, fr) {
    ## based on mkBlist

    noGrpFac <- is.null(lme4::findbars(bar))

    if(!noGrpFac) {
        linFormLang <- bar[[2]] # language object specifying linear model
            grpLang <- bar[[3]] # language object specifying grouping factor

        fr <- lme4::factorize(bar, fr)
        nm <- deparse(grpLang)
        ## try to evaluate grouping factor within model frame ...
        if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac),
                                                   list(fac = grpLang)), fr),
                                   error = function(e) NULL)))
            stop("couldn't evaluate grouping factor ",
                 nm, " within model frame:",
                 " try adding grouping factor to data ",
                 "frame explicitly if possible", call. = FALSE)
        if (all(is.na(ff)))
            stop("Invalid grouping factor specification, ",
                 nm, call. = FALSE)
    } else { # noGrpFac
        linFormLang <- bar
        ff <- nm <- NA
    }
    mm <- model.matrix(eval(substitute( ~ foo, list(foo = linFormLang))), fr)
    return(list(modMat = mm, grpFac = ff, grpName = nm))
}

##' Simplify factor list over random effects terms
##'
##' @param facList list of grouping factors over random effects terms
##' @return collapse repeated factors and add an \code{assign}
##' attribute for indicating which factors relate to which terms
##' @export
simplifyFacList <- function(facList) {
    fnms <- names(facList)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        facList <- facList[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    } else asgn <- seq_along(facList)
    names(facList) <- ufn
    facList <- do.call(data.frame, c(facList, check.names = FALSE))
    attr(facList, "assign") <- asgn
    return(facList)
}

##' Generate phylogenetic test data
##'
##' @param seed random seed
##' @param n number of sites
##' @param m number of species
##' @param form formula with \code{y} (comm dat), \code{x} (env),
##' \code{z} (trait), \code{species}, or \code{sites}
##' @param power power for \code{Grafen} method
##' @param covarSim covariance parameters
##' @param fixefSim fixed effect parameters
##' @param perm permute data list?
##' @importMethodsFrom Matrix t
##' @export
simTestPhyloDat <- function(seed = 1, n = 10, m = 30,
                            form = y ~ 1 + (1 | species),
                            power = 0.1,
                            covarSim = 1, fixefSim,
                            perm = FALSE) {
    set.seed(seed)
    dl <- dims_to_vars(data.list(y = 1 * (matrix(rnorm(n * m), n, m) > 0),
                                 x = rnorm(n), z = rnorm(m),
                                 dimids = c("sites", "species")))
    if(perm) dl <- aperm(dl, c(2, 1))
    df <- as.data.frame(dl)
    phy <- ape::rtree(n = m)
    phy <- ape::compute.brlen(phy, method = "Grafen", power = power)
    Vphy <- stanCov(ape::vcv(phy))
    dimnames(Vphy) <- rep(list(1:m), 2)
    covList <- list(species = Vphy)
    parsedForm <- glmercFormula(form, df, covList = covList, strList = list())
    parsedForm <- within(parsedForm, Lambdat@x[] <- mapToCovFact(covarSim))
    X <- model.matrix(nobars(form), df) # fixed effects design matrix
    Z <- t(parsedForm$Lambdat %*% parsedForm$Zt) # random effects design
                                        # matrix with
                                        # phylogenetic
                                        # covariances
    u <- rnorm(ncol(Z)) # whitened random effects
    if(missing(fixefSim)) fixefSim <- rnorm(ncol(X))
    p <- plogis(as.numeric(X %*% fixefSim + Z %*% u)) # probability of observation
    dl$y <- rbinom(nrow(df), 1, p) # presence-absence data
    dimnames(dl)[[2]] <- phy$tip.label
    return(list(dl = dims_to_vars(dl),
                cv = Vphy,
                ph = phy))
}


##' Inverse of the n-choose-2 function
##'
##' @param m a vector coercible to integer
##' @export
##' @examples
##' nChoose2Inv(choose(1:10, 2)) # = 1:10
nChoose2Inv <- function(m) as.integer((sqrt(1 + 8 * as.integer(m)) + 1)/2)




## Modified from lme4's modification of stats simulate() functions

gaussian_simfun <- function(wts, nsim, ftd) {
    if (any(wts != 1)) warning("ignoring prior weights")
    rnorm(nsim*length(ftd), ftd, sd=sigma(object))
}

binomial_simfun <- function(wts, nsim, ftd=fitted(object)) {
    rbinom(nsim, size = wts, prob = ftd)/wts
}

poisson_simfun <- function(wts, nsim, ftd) {
        ## A Poisson GLM has dispersion fixed at 1, so prior weights
        ## do not have a simple unambiguous interpretation:
        ## they might be frequency weights or indicate averages.
        if (any(wts != 1)) warning("ignoring prior weights")
        rpois(nsim*length(ftd), ftd)
    }


## FIXME: need a gamma.shape.merMod method in order for this to work.
##        (see initial shot at gamma.shape.merMod below)
Gamma_simfun <- function(object, nsim, ftd=fitted(object)) {
    stop("not implemented")
    wts <- weights(object)
    if (any(wts != 1)) message("using weights as shape parameters")
    ## ftd <- fitted(object)
    shape <- MASS::gamma.shape(object)$alpha * wts
    rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
}

gamma.shape.merMod <- function(object, ...) {
    stop("not implemented")
    if(family(object)$family != "Gamma")
	stop("Can not fit gamma shape parameter because Gamma family not used")

    y <- getME(object, "y")
    mu <- getME(object, "mu")
    w <- weights(object)
                                        # Sec 8.3.2 (MN)
    L <- w*(log(y/mu)-((y-mu)/mu))
    dev <- -2*sum(L)
                                        # Eqs. between 8.2 & 8.3 (MN)
    Dbar <- dev/length(y)
    structure(list(alpha = (6+2*Dbar)/(Dbar*(6+Dbar)),
		   SE = NA), # FIXME: obtain standard error
	      class = "gamma.shape")
}


## FIXME: include without inducing SuppDists dependency?
## inverse.gaussian_simfun <- function(object, nsim, ftd=fitted(object)) {
##     if(is.null(tryCatch(loadNamespace("SuppDists"),
##                         error = function(e) NULL)))
##         stop("need CRAN package 'SuppDists' for the 'inverse.gaussian' family")
##     wts <- weights(object)
##     if (any(wts != 1)) message("using weights as inverse variances")
##     SuppDists::rinvGauss(nsim * length(ftd), nu = ftd,
##                          lambda = wts/summary(object)$dispersion)
## }

## in the original MASS version, .Theta is assigned into the environment
## (triggers a NOTE in R CMD check)
negative.binomial_simfun <- function (object, nsim, ftd=fitted(object))
{
    stop("not implemented yet")
    ## val <- rnbinom(nsim * length(ftd), mu=ftd, size=.Theta)
}

simfunList <- list(gaussian = gaussian_simfun,
		   binomial = binomial_simfun,
		   poisson  = poisson_simfun,
		   Gamma    = Gamma_simfun,
		   negative.binomial = negative.binomial_simfun)


##' Convert row and column indices to vector indices, for a triangular
##' matrix
##'
##' @param rowInds,colInds vectors of 1-based row and column indices
##' @param maxInd maximum index (could be larger than
##' \code{max(rowInds, colInds)})
##' @param type type of parameterization (currently only
##' \code{"distClass"} available)
##' @export
##' @examples
##' n <- 5
##' set.seed(1)
##' X <- dist(matrix(rnorm((n + 1) * 2), n + 1, 2))
##' rowInds <- rep(1:n, 1:n)
##' colInds <- sequence(1:n)
##' X[triInds(rowInds, colInds, n)]
##' X
triInds <- function(rowInds, colInds, maxInd,
                    type = c("distClass")) {
    stopifnot(type == "distClass")
    rowInds + (colInds - 1) * maxInd - choose(colInds, 2)
}


##' Count unique values
##'
##' @param x numeric (for \code{countUnique}) or integer (for
##' \code{countInRange}) vector
##' @return data frame with \code{uniqueVals} and \code{counts}
##' columns
##' @rdname count
##' @export
countUnique <- function(x) {
    sx <- sort(x)
    counts <- diff(c(which(!duplicated(sx)), length(x) + 1))
    data.frame(uniqueVals = unique(sx),
               counts = counts)
}

##' @rdname count
##' @export
countInRange <- function(x) {
    x <- as.integer(x)
    vals <- min(x):max(x)
    counts <- sapply(vals, function(xx) sum(xx == x))
    data.frame(vals = vals, counts = counts)
}

##' Flatten an integer vector
##'
##' @param x vector coercible to nonnegative integers
##' @export
flattenIntVec <- function(x) {
    x <- abs(as.integer(x))
    match(x, sort(unique(x)))
}

##' Assign to an environment the value of an expression evaluated in
##' another environment (or list or data frame)
##'
##' @param expr expression to evaluate
##' @param name name of object in \code{assignEnv}
##' @param data environment (or list or data frame) in which to
##' evaluate \code{expr}
##' @param envir environment in which to store the results of
##' \code{expr} as \code{name}
##' @param enclos enclosing environemtn to be passed to
##' \code{\link{assign}}
##' @export
assignWith <- function(expr, name, data, envir, enclos = parent.frame()) {
    assign(name, eval(substitute(expr), data, enclos = enclos), envir = envir)
}


mapplyInvList <- function(vecList, lens) {
    matList <- lapply(mapply(matrix, vecList, lens, lens, SIMPLIFY = FALSE), t)
    diagList <- mapply(diag, nrow = lens, ncol = lens, SIMPLIFY = FALSE)
    mapply(backsolve, matList, diagList, SIMPLIFY = FALSE)
}



##' Show skeleton for some type of function
##'
##' @param whichSkel length-one character giving the name of the
##' skeleton
##' @param package the default is probably what you want
##' @export
showSkeleton <- function(whichSkel, package = "lme4ord") {
   fn <- system.file("skeletons",
                      paste(whichSkel, "R", sep = "."),
                      package = package)
   if(nchar(fn) == 0L) {
       availableSkels <- list.files(system.file("skeletons", package = "lme4ord"))
       nc <- nchar(availableSkels)
       avsk <- sapply(strsplit(availableSkels, "\\."),
                      function(xx) paste(xx[-length(xx)], collapse = ""))
       stop("skeleton not found. ",
            "here are the available skeletons: \n",
            avsk)
   }
   cat(paste(c(readLines(fn), ""), collapse = "\n"))
}


.safeExtractDiag <- function(x) {
    if(length(x) == 1) return(x)
    diag(x)
}

##' Transpose a list
##'
##' @param lst \code{\link{list}}
##' @export
listTranspose <- function(lst) {
    lstExtract <- function(i) lapply(lst, "[[", i)
    setNames(lapply(names(lst[[1]]), lstExtract), names(lst[[1]]))
}


subRagByLens <- function(x, lens) {
    split(x, rep(seq_along(lens), lens)) ## split no good ... order of levels !
}

##' Indices for square dense matrices
##'
##' @param n matrix size
##' @rdname denseInds
##' @export
denseDiagInds <- function(n) {
    seq(1, n^2, length.out = n)
}

##' @rdname denseInds
##' @export
denseLowerInds <- function(n) {
    unlist(lapply(1:(n-1), function(i) (i - 1) * n + ((i + 1):n)))
}

##' @rdname denseInds
##' @export
denseUpperInds <- function(n) {
    unlist(lapply(1:(n-1), function(i) i * n + (1:i)))
}



##----------------------------------------------------------------------
## edgeNodeTools
##----------------------------------------------------------------------

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

##' @importFrom ape as.phylo
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

##' Orthogonal procrustean rotation matrix
##' 
##' @param X input matrix
##' @param Z target matrix
##' @return the rotation matrix
##' @export
orthProcrustesRotMat <- function(X, Z) {
    sol <- svd(crossprod(Z, X))
    return(sol$v %*% t(sol$u))
}
