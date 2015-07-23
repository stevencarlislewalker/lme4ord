## ----------------------------------------------------------------------
## Definition of an lme4ord utility:
##
## "A function that does _not_ make use of functions from any other
## file in the package."
##
## The purpose of this definition is to preserve modularity.
## ----------------------------------------------------------------------

## ----------------------------------------------------------------------
## strucGlmer object extraction functions
## ----------------------------------------------------------------------

##' Make parameter indices
##'
##' One faces two somewhat opposing design goals when deciding on how
##' to represent the parameters of structured generalized linear mixed
##' models. On one hand, it would be nice to organize the
##' representation in terms of the four different types of parameters
##' (see below): (\code{covar}, \code{fixef}, \code{loads},
##' \code{weigh}). On the other hand, the nonlinear optimizer takes
##' vector valued parameter sets. The \code{mkParInds} function is
##' used to get the best of both worlds, by constructing a list of
##' indices for extracting various types of parameters from a
##' parameter vector.
##' 
##' The \code{\link{lme4ord}} package keeps track of parameters using
##' a simple named list of parameter vectors. The names correspond to
##' different types of parameters, described in the following list:
##' \describe{
##' 
##' \item{\code{covar}}{Parameters determining the transposed relative
##' covariance factor, \code{Lambdat}.}
##'
##' \item{\code{fixef}}{The fixed effects coefficients.}
##'
##' \item{\code{loads}}{Parameters (called loadings) determining the
##' random effects model matrix, \code{Zt}.}
##'
##' \item{\code{weigh}}{Parameters determining the observation weights
##' (experimental).}
##'
##' }
##' \code{mkParInds} constructs the indices required to extract the
##' different types of parameters in \code{unlist(parList)}.
##'
##' @return A list with the same names as \code{parList} with the
##' indices for each type of parameter.
##' 
##' @param parList named list of parameters with possible names (see
##' details): (\code{covar}, \code{fixef}, \code{loads}, \code{weigh})
##' @export
##' @examples
##' set.seed(1)
##' parList <- list(covar = 1, fixef = c(0, 0), loads = rnorm(10))
##' parInds <- mkParInds(parList)
##' parVec <- unlist(parList)
##' parVec[parInds$fixef]
mkParInds <- function(parList) {
    if(!is.recursive(parList)) stop("parList must be a list")
    if(length(parList) == 1L) return(lapply(parList, seq_along))
    parInds <- mapply(`+`,
                      lapply(parList, seq_along),
                      c(0, cumsum(lapply(parList, length))[-length(parList)]),
                      SIMPLIFY = FALSE)
    names(parInds) <- names(parList) ## too paranoid?
    keepers <- sapply(parInds, length) > 0
    parInds[keepers]
}

## ----------------------------------------------------------------------
## Initial values -- get and set init parameter vectors for repeated
## sparse matrices
## ----------------------------------------------------------------------

##' Get and set initial parameter values for repeated sparse matrices
##'
##' The initial values for a repeated sparse matrix are stored in the
##' environment of its transformation function, where they have the
##' name \code{init}.  Random effects term structures
##' (\code{\link{reTrmStruct}}) contain two repeated sparse matrices:
##' the transposed relative covariance factor, \code{Lambdat}, and the
##' transposed random effects model matrix, \code{Zt}. There are
##' \code{getReTrm.reTrmStruct} and \code{setReTrm.reTrmStruct}
##' methods for getting and setting the initial values of these
##' matrices directly.  Setting initial values does not change the
##' values themselves, only the initial values. Use the
##' \code{\link{update.repSparse}} methods to actually update the
##' values to the initial values.
##' 
##' @param x object
##' @param ... not yet used
##' @rdname getInit
##' @export
##' @examples
##' set.seed(1)
##' m1 <- as.repSparse(matrix(rnorm(6), 2, 3))
##' m2 <- repSparseCompSymm(1.2, -0.2, 5)
##' getInit(m1)
##' getInit(m2)
##' setInit(m2, c(10, 0))
getInit <- function(x, ...) UseMethod("getInit")

##' @rdname getInit
##' @export
getInit.default <- function(x, ...) x$init

##' @rdname getInit
##' @export
getInit.repSparse <- function(x, ...) getInit(x$trans)

##' @rdname getInit
##' @export
getInit.function <- function(x, ...) environment(x)$init

##' @rdname getInit
##' @export
getInit.reTrmStruct <- function(x, ...) {
    list(initCovar = getInit(x$Lambdat),
         initLoads = getInit(x$Zt))
}

##' @rdname getInit
##' @export
setInit <- function(x, ...) UseMethod("setInit")

##' @param init initial value for the parameter vector
##' @rdname getInit
##' @export
setInit.default <- function(x, init, ...) {
    assign("init", init, envir = as.environment(x))
}

##' @rdname getInit
##' @export
setInit.repSparse <- function(x, init, ...) {
    assign("init", init, envir = environment(x$trans))
}

##' @rdname getInit
##' @export
setInit.function <- function(x, init, ...) {
    assign("init", init, envir = environment(x))
}

##' @param initCovar initial value for the covariance parameter vector
##' @param initLoads initial value for the loadings parameter vector
##' @rdname getInit
##' @export
setInit.reTrmStruct <- function(x, initCovar, initLoads, ...) {
    setInit(x$Lambdat, initCovar)
    setInit(x$Zt,      initLoads)
}

##' @param ll length of parameter vector
##' @rdname getInit
##' @export
parLength <- function(ll) {
    initll <- getInit(ll)
    if(is.null(initll) || all(is.na(initll))) return(0L)
    return(length(initll))
}

##' Extract factors
##'
##' Intended for factor analysis (\code{\link{factAnal}}), and
##' modelled after the \code{\link{loadings}} function in the
##' \code{stats} package.
##'
##' @param x an object with \code{x$factors}
##' @param ... not used
##' @export
factors <- function(x, ...) x$factors

## ----------------------------------------------------------------------
## Phylogenetic 

##' Standardize covariance matrix to determinant one
##'
##' @param covMat covariance matrix
##' @export
stanCov <- function(covMat) {
    covMat / (det(covMat)^(1/nrow(covMat)))
}


##' Inverse of the n-choose-2 function
##'
##' @note No checks are made for non-integer input or output.
##'
##' @param m a vector coercible to integer
##' @export
##' @examples
##' nChoose2Inv(choose(1:10, 2)) # = 1:10
nChoose2Inv <- function(m) (sqrt(1 + 8 * m) + 1)/2


## ----------------------------------------------------------------------
## Modified from lme4's modification of stats simulate() functions

##' Family simulation functions
##'
##' @param object typically a \code{\link{family}} object
##' @param ... additional objects to pass on (not yet used).
##' @return a simulation function taking three arguments: \code{wts},
##' \code{nsim}, and \code{ftd}.
##' @export
##' @examples
##' set.seed(1)
##' simFun <- familySimFun(binomial)
##' simFun(wts = rep(1, 10),
##'        nsim = 10,
##'        ftd = binomial()$linkinv(rnorm(10)))
familySimFun <- function(object, ...) {
    UseMethod("familySimFun")
}

##' @rdname familySimFun
##' @export
familySimFun.family <- function(object, ...) {
    dummyClass <- paste(object$family, "Family", sep = "")
    familySimFun(structure(list(), class = dummyClass))
}

##' @rdname familySimFun
##' @export
familySimFun.character <- function(object, ...) {
    dummyClass <- paste(object, "Family", sep = "")
    familySimFun(structure(list(), class = dummyClass))
}

##' @rdname familySimFun
##' @export
familySimFun.function <- function(object, ...) familySimFun(object())

##' @rdname familySimFun
##' @export
familySimFun.gaussianFamily <- function(object, ...) {
    function(wts, nsim, ftd) {
        if (any(wts != 1)) warning("ignoring prior weights")
        rnorm(nsim*length(ftd), ftd, sd = 1) ## FIXME: sd = sigma(object)
    }
}

##' @rdname familySimFun
##' @export
familySimFun.binomialFamily <- function(object, ...) {
    function(wts, nsim, ftd=fitted(object)) {
        rbinom(nsim, size = wts, prob = ftd)/wts
    }
}

##' @rdname familySimFun
##' @export
familySimFun.poissonFamily <- function(object, ...) {
    function(wts, nsim, ftd) {
        ## A Poisson GLM has dispersion fixed at 1, so prior weights
        ## do not have a simple unambiguous interpretation:
        ## they might be frequency weights or indicate averages.
        if (any(wts != 1)) warning("ignoring prior weights")
        rpois(nsim*length(ftd), ftd)
    }
}

## FIXME: need a gamma.shape.merMod method in order for this to work.
##        (see initial shot at gamma.shape.merMod below)

##' @rdname familySimFun
##' @export
familySimFun.GammaFamily <- function(object, ...) {
    function(object, nsim, ftd=fitted(object)) {
        stop("not implemented")
        wts <- weights(object)
        if (any(wts != 1)) message("using weights as shape parameters")
        ## ftd <- fitted(object)
        shape <- MASS::gamma.shape(object)$alpha * wts
        rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
    }
}

gamma.shape.merMod <- function(object, ...) {
    stop("not implemented")
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
## negative.binomial_simfun <- function (object, nsim, ftd=fitted(object))
## {
##     stop("not implemented yet")
##     ## val <- rnbinom(nsim * length(ftd), mu=ftd, size=.Theta)
## }

## simfunList <- list(gaussian = gaussian_simfun,
## 		   binomial = binomial_simfun,
## 		   poisson  = poisson_simfun,
## 		   Gamma    = Gamma_simfun,
## 		   negative.binomial = negative.binomial_simfun)


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
##' @return \code{\link{data.frame}} with \code{uniqueVals} and
##' \code{counts} columns
##' @rdname count
##' @export
##' @examples
##' x <- rep(1:5, 5:1)
##' countUnique(x)
##' countInRange(x)
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
##' Map integers in a vector, \code{x}, to the integers from \code{1,
##' ..., n}, where \code{n} is the number of unique integers in
##' \code{x}.
##'
##' @param x vector coercible to nonnegative integers
##' @export
##' @examples
##' flattenIntVec(rep(seq(0, 100, length = 5), 1:5))
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
##' @examples
##' showSkeleton("setReTrm")
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
##' @examples
##' set.seed(1)
##' X <- repSparse(c(1, 2, 1, 2), c(1, 1, 2, 2), 1:4, rnorm(4))
##' Y <- repSparse(c(1, 2, 1), c(1, 1, 2), 1:3, rnorm(3))
##' Z <- repSparse(c(1, 2), c(1, 2), 1:2, rnorm(2))
##' listTranspose(list(X = X, Y = Y, Z = Z))
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
##' @examples
##' set.seed(1)
##' A <- matrix(rnorm(25), 5, 5)
##' A[denseDiagInds(5)]  ## = diag(A)
##' A[denseLowerInds(5)] ## = A[lower.tri(A)]
##' A[denseUpperInds(5)] ## = A[upper.tri(A)]
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
##' if (require("ape")) {
##'    phy <- rtree(5)
##'    plot(phy)
##'    edgelabels()
##'    edgeTipIndicator(phy)
##'    (ee <- edgeTipIndicator(phy))
##'    if (require("Matrix")) {
##'        image(Matrix(ee),sub="",xlab="tips",ylab="branches")
##'    }
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

##' Exponentially decaying covariances
##'
##' Compute the exponential decay associated with an
##' \code{\link{expDecay}} term.
##'
##' @param matPars parameters of an \code{\link{expDecay}} term
##' (FIXME: give more info).
##' @param minCov,distCutoff see \code{\link{expDecay}}.
##' @param nPoints number of points at which to evaluate the curve.
##' @export
##' @examples
##' covExpDecay(c(1, 1))
covExpDecay <- function(matPars, minCov = 1e-3, distCutoff = 2, nPoints = 100) {
    edgeDists <- seq(0, distCutoff, length = nPoints)
    q1 <- (minCov - 1)/(exp(-(matPars[1]) * distCutoff) - 1)
    q2 <- 1 - q1
    edgeCovs <- matPars[2] * (q2 + q1 * exp(-(matPars[1]) * edgeDists))
    return(data.frame(edgeDists = edgeDists, edgeCovs = edgeCovs))
}

##----------------------------------------------------------------------
## factor analysis tools
##----------------------------------------------------------------------

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

## ----------------------------------------------------------------------
## formula parsing tools
## ----------------------------------------------------------------------

##' Simplify factor list over random effects terms (not currently
##' used)
##'
##' @param facList list of grouping factors over random effects terms
##' @return collapse repeated factors and add an \code{assign}
##' attribute for indicating which factors relate to which terms
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


## ----------------------------------------------------------------------
## print utilities copied from lme4 (FIXME: export from lme4???)
## ----------------------------------------------------------------------

.prt.methTit <- function(mtit, class) {
    if(nchar(mtit) + 5 + nchar(class) > (w <- getOption("width"))) {
	## wrap around
	mtit <- strwrap(mtit, width = w - 2, exdent = 2)
	cat(mtit, " [",class,"]", sep = "", fill = TRUE)
    } else ## previous: simple one-liner
	cat(sprintf("%s ['%s']\n", mtit, class))
}

.prt.family <- function(famL) {
    if (!is.null(f <- famL$family)) {
	cat.f(" Family:", f,
	      if(!is.null(ll <- famL$link)) paste(" (", ll, ")"))
    }
}

.prt.resids <- function(resids, digits, title = "Scaled residuals:", ...) {
    cat(title,"\n")
    ## FIXME: need testing code
    rq <- setNames(zapsmall(quantile(resids, na.rm=TRUE), digits + 1L),
                   c("Min", "1Q", "Median", "3Q", "Max"))
    print(rq, digits = digits, ...)
    cat("\n")
}

.prt.call <- function(call, long = TRUE) {
    if (!is.null(cc <- call$formula))
	cat.f("Formula:", deparse(cc))
    if (!is.null(cc <- call$data))
	cat.f("   Data:", deparse(cc))
    if (!is.null(cc <- call$weights))
        cat.f("Weights:", deparse(cc))
    if (!is.null(cc <- call$offset))
        cat.f(" Offset:", deparse(cc))
    if (long && length(cc <- call$control) &&
	!identical((dc <- deparse(cc)), "lmerControl()"))
	## && !identical(eval(cc), lmerControl()))
	cat.f("Control:", dc)
    if (!is.null(cc <- call$subset))
	cat.f(" Subset:", deparse(cc))
    }

getLlikAIC <- function(object, cmp = object@devcomp$cmp) {
    llik <- logLik(object)   # returns NA for a REML fit - maybe change?
    AICstats <- {
	if(isREML(object)) cmp["REML"] # *no* likelihood stats here
	else {
	    c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
	      deviance = lme4:::devCritFun(object),
              df.resid = df.residual(object))
	}
    }
    list(logLik = llik, AICtab = AICstats)
}

.prt.aictab <- function(aictab, digits = 1) {
    t.4 <- round(aictab, digits)
    if (length(aictab) == 1 && names(aictab) == "REML")
	cat.f("REML criterion at convergence:", t.4)
    else {
        ## slight hack to get residual df formatted as an integer
        t.4F <- format(t.4)
        t.4F["df.resid"] <- format(t.4["df.resid"])
        print(t.4F, quote = FALSE)
    }
}

.prt.VC <- function(varcor, digits, comp, formatter = format, ...) {
    cat("Random effects:\n")
    fVC <- if(missing(comp))
	lme4:::formatVC(varcor, digits = digits, formatter = formatter)
    else
	lme4:::formatVC(varcor, digits = digits, formatter = formatter, comp = comp)
    print(fVC, quote = FALSE, digits = digits, ...)
}

.prt.grps <- function(ngrps, nobs) {
    cat(sprintf("Number of obs: %d, groups: ", nobs),
        paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "),
        fill = TRUE)
}

## FIXME: print header ("Warnings:\n") ?
##  change position in output? comes at the very end, could get lost ...
.prt.warn <- function(optinfo, summary=FALSE, ...) {
    ## check all warning slots: print numbers of warnings (if any)
    cc <- optinfo$conv$opt
    msgs <- unlist(optinfo$conv$lme4$messages)
    ## can't put nmsgs/nwarnings compactly into || expression
    ##   because of short-circuiting
    nmsgs <- length(msgs)
    warnings <- optinfo$warnings
    nwarnings <- length(warnings)
    if (cc>0 || nmsgs>0 || nwarnings>0) {
        if (summary) {
            cat(sprintf("convergence code %d; %d optimizer warnings; %d lme4 warnings",
                cc,nmsgs,nwarnings),"\n")
        } else {
            cat(sprintf("convergence code: %d",cc),
                msgs,
                unlist(warnings),
                sep="\n")
            cat("\n")
        }
    }
}

cat.f <- function(...) cat(..., fill = TRUE)

famlink <- function(object, resp = object@resp) {
    if(is(resp, "glmResp"))
	resp$family[c("family", "link")]
    else list(family = NULL, link = NULL)
}

