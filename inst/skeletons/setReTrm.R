##' @param object an object of class MYCLASS constructed by mkReTrmStructs
##' @param addArgsList a list of objects passed to strucGlmer, which
##' might be required to process any additional arguments passed
##' through structured random effects terms in the formula passed to
##' strucGlmer
##' @param devfunEnv the environment of the deviance function
##' (probably best to ignore this at first). used to access model
##' objects when updating the parameters associated with this term
##' (e.g. when the covariance factor depends on the fitted values)
setReTrm.MYCLASS <- function(object, addArgsList,
                             devfunEnv = NULL) {

    ## ----------------------------------------------------------------------
    ## Additional arguments
    ## 
    ## uncomment to get any additional arguments passed through the
    ## formula
    ##
    ## note: object$addArgs are the additional arguments passed
    ## through the formula
    ## 
    ## note: the -1 simply removes the class name, which shouldn't be
    ## an argument.
    ## ----------------------------------------------------------------------
    
    ## addArgs <- getAddArgs(object$addArgs[-1], addArgsList)


    ## ----------------------------------------------------------------------
    ## Transposed model matrix (as an object of class repSparse)
    ##
    ## below gives the strategy for a standard lme4 term, which is a
    ## khatri-rao product between the raw model matrix and an
    ## indicator matrix for the levels of the grouping factor
    ##
    ## note: also can be used as a loadings matrix for factor analytic
    ## models
    ##
    ## note: unless a factor analytic term is desired, it is necessary
    ## to wrap the matrix in a resetTransConst so that the matrix is
    ## not parameterized
    ##
    ## note: row binded with the transposed model matrix for the other
    ## terms to obtain the full transposed random effects model matrix
    ## ----------------------------------------------------------------------
    
    Zt <- resetTransConst(kr(t(as.repSparse(object$modMat)),
                               as.repSparse(object$grpFac)))
    

    ## ----------------------------------------------------------------------
    ## Relative covariance factor (as an object of class repSparse)
    ##
    ## below gives the strategy for a standard lme4 term, which is a
    ## block-diagonal lower triangular matrix, with one block per
    ## level of the grouping factor
    ##
    ## note: block-diagonal binded with the relative covariance
    ## factors for the other terms to obtain the full model relative
    ## covariance factor
    ## ----------------------------------------------------------------------

    nc <-    ncol(object$modMat)
    nl <- nlevels(object$grpFac)
    templateBlock <- repSparseTri(   diagVals = rep(1,        nc    ),
                                  offDiagVals = rep(0, choose(nc, 2)),
                                  low = FALSE)
    Lambdat <- rep(templateBlock, nl, type = "diag")

    ## ----------------------------------------------------------------------
    ## Set optimizer bounds (optional)
    ##
    ## below is a silly example setting bounds between -1 and 1 for
    ## all parameters
    ##
    ## note: by default standard lme4 lower and upper bounds are used
    ## which are: ifelse(init, 0, -Inf) and rep(Inf, length(init))
    ## ----------------------------------------------------------------------
    
    ## lowerLoads <- rep(-1, length(getInit(Zt)))
    ## lowerCovar <- rep(-1, length(getInit(Lambdat)))
    ## upperLoads <- rep( 1, length(getInit(Zt)))
    ## upperCovar <- rep( 1, length(getInit(Lambdat)))

    ## ----------------------------------------------------------------------
    ## Package up the results
    ##
    ## note: can pass any specialized bounds here (see ?packReTrm)
    ## ----------------------------------------------------------------------

    packReTrm(object, Zt, Lambdat)
}

