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
    
    fr <- lme4:::factorize(bar, fr)
    grpLang <- bar[[3]]
    linFormLang <- bar[[2]]
    nm <- deparse(grpLang)
    ## try to evaluate grouping factor within model frame ...
    if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac),
                                               list(fac = grpLang)), fr),
                error=function(e) NULL)))
        stop("couldn't evaluate grouping factor ",
             nm," within model frame:",
             " try adding grouping factor to data ",
             "frame explicitly if possible",call.=FALSE)
    if (all(is.na(ff)))
        stop("Invalid grouping factor specification, ",
             nm, call. = FALSE)
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
    phy <- ape:::rtree(n = m)
    phy <- ape:::compute.brlen(phy, method = "Grafen", power = power)
    Vphy <- stanCov(ape:::vcv(phy))
    dimnames(Vphy) <- rep(list(1:m), 2)
    covList <- list(species = Vphy)
    parsedForm <- glmercFormula(form, df, covList = covList)
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
