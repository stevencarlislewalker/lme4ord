glmerfFormula <- function(formula, data = NULL, family = binomial,
                          latentName = "latent", latentDims = 0L,
                          loadingsDim = 1L, loadingPen = 0L,
                          loadingsPat = NULL, ...) {
    
}

##' Modular functions for generalized bilinear mixed models
##' @export
mkGblmerDevfun <- function(fr, X, linReTrms, bilinReTrms,
                           family, nAGQ = 1, verbose = 0L, ...) {
    reTrms <- joinReTrms(bilinReTrms, linReTrms)
    dfun0  <- mkGlmerDevfun(fr, X, reTrms, family(), nAGQ, verbose, ...)
    dfun   <- updateGlmerDevfun(dfun0, reTrms)
    rho    <- environment(dfun)

    rho$q     <- reTrms$q
    rho$nth   <- reTrms$nth
    rho$nCnms <- reTrms$nCnms
                                        # which Zt@x elements
                                        # represent loadings
    rho$Zwhich <- rho$pp$Zt@i %in% (seq_len(rho$q[1]) - 1)
                                        # mapping from loadings to the
                                        # Zt@x elements that represent
                                        # loadings
    rho$Zind <- with(rho, pp$Zt@x[Zwhich])
    rho$Ztx <- rho$pp$Zt@x
    
    return(dfun)
}
