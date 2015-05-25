##' \code{lme4ord}: community ecological mixed effects models with \code{lme4}
##'
##' The idea behind lme4ord is to provide a flexible platform (and
##' maybe one day an API) for using lme4 computational machinery to
##' fit generalized linear mixed models to ecological community
##' data. The reason for this package is that I have found that (1)
##' community data requires more flexibility than the standard lme4
##' interface can provide but that (2) lme4 can usually be hacked to
##' fit the required model. Although I am a community ecologist, it is
##' likely that lme4ord is more broadly useful. However, as a
##' community ecologist am mostly interested in the following types of
##' models:
##' \describe{
##' \item{1}{Factor analytic and structural equation latent variable models}
##' \item{2}{Models with phylogenetic correlations}
##' \item{3}{Models with spatial correlations}
##' }
##' @docType package
##' @name lme4ord
##' @aliases lme4ord package-lme4ord lme4ord-package
##' @import MASS lattice minqa
NULL
