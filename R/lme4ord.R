##' \code{lme4ord}: Generalized Linear Mixed Models with Flexible
##' (Co)Variance Structures
##'
##' Generalized linear mixed models with flexible (co)variance
##' structures. The package was written to solve applied problems in
##' biology and ecology, but should be more general than that.
##'
##' @section Motivation:
##' The original goal of lme4ord is to provide a flexible platform and
##' API, for using the computational machinery from the \code{lme4}
##' package to fit generalized linear mixed models to ecological
##' community data. Models for community data typically require more
##' flexibility than the standard \code{lme4} interface can provide,
##' despite the fact that \code{lme4} can usually be hacked to fit
##' such models. \code{lme4ord} was developed to avoid such hacking by
##' providing infrastructure for specifying a wide range of
##' models. For example, \code{lme4ord} is able to fit the following
##' types of models:
##' \itemize{
##'   \item Factor analytic and structural equation latent variable models
##'   \item Models with phylogenetic correlations
##'   \item Models with spatial correlations
##' }
##' Despite this origin in community ecology, \code{lme4ord} is more
##' broadly useful for fitting a very wide variety of structured
##' generlized linear mixed models.
##' @docType package
##' @name lme4ord
##' @aliases lme4ord package-lme4ord lme4ord-package
##' @import MASS lattice minqa
NULL
