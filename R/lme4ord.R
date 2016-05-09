##' \code{lme4ord}: Generalized Linear Mixed Models with Flexible
##' (Co)Variance Structures
##'
##' @section Motivation:
##' 
##' The \code{lme4ord} package was written to solve applied problems
##' in biology and ecology, but should be more general than that.  The
##' original goal of \code{lme4ord} is to provide a flexible platform
##' and API, for using the computational machinery from the
##' \code{lme4} package to fit generalized linear mixed models to
##' ecological community data. Models for community data typically
##' require more flexibility than the standard \code{lme4} interface
##' can provide, despite the fact that \code{lme4} can usually be
##' hacked to fit such models. \code{lme4ord} was developed to avoid
##' such hacking by providing infrastructure for specifying a wide
##' range of models. For example, \code{lme4ord} is able to fit the
##' following types of models:
##' 
##' \itemize{
##' 
##' \item Factor analytic and structural equation latent variable
##' models (\code{\link{factAnal}})
##'
##' \item Models with phylogenetic correlations
##' (\code{\link{phyloEdge}})
##'
##' \item Models with spatial correlations (\code{\link{expDecay}})
##'
##' \item Any \code{\link{corStruct}} model from the \code{nlme}
##' package (\code{\link{nlmeCorStruct}})
##'
##' }
##' 
##' Despite this origin in community ecology, \code{lme4ord} is more
##' broadly useful for fitting a very wide variety of structured
##' generlized linear mixed models.
##'
##' @section Modularization:
##'
##' The design of \code{lme4ord} is intended to facilitate use by
##' people with a range of skill levels and interests. In order to
##' achieve this goal, \code{lme4ord} is composed of several modules
##' (illustrated in the figure below). Those who choose to interact
##' with all modules will benefit from maximum flexibility, whereas
##' more simple uses may require only a few modules. 
##'
##' \figure{modularizationDiagram.png}
##'
##' Arrows point from modules that make use of functions from the
##' module to which they point.
##'
##' \describe{
##'
##' \item{structured glmm}{The model fitting module consists of a
##' single function, \code{\link{strucGlmer}}. This function calls the
##' next four modules in the order given.}
##'
##' \item{formula parsing}{The formula parsing module takes a mixed
##' model formula and data, and returns a set of objects required to
##' fit the specified model. This module may be accessed directly via
##' the \code{\link{strucParseFormula}} function. In addition to
##' standard \code{lme4}-style random effects terms, mixed model
##' formulas may consist of special \code{\link{reTrmStruct}} object
##' constructed using the random-effects structures module described
##' below. The \code{\link{findReTrmClasses}} function may be used to
##' list all such special structures that are available. Also see:
##' \code{\link{splitForm}};
##' \code{\link{simulate.strucParseFormula}}.}
##'
##' \item{deviance function construction}{This module is what links a
##' \code{\link{strucParseFormula}} object with the underlying
##' \code{C++} machinery of \code{lme4}. The function
##' \code{\link{mkGeneralGlmerDevfun}} may be used to construct a
##' function that uses \code{C++} code to compute the deviance (minus
##' twice the log-likelihood) for a given set of parameters. Also see
##' \code{\link{mkPenaltyFun}}.}
##'
##' \item{optimization}{Currently \code{lme4ord} uses the
##' \code{\link{bobyqa}} function in the \code{minqa} package to
##' optimize deviance functions. }
##'
##' \item{output}{There are a variety of methods and extractor
##' functions for exploring fitted models of
##' \code{\link{strucGlmer-class}}. These objects are constructed by
##' \code{\link{strucGlmer}}, which calls
##' \code{\link{mkStrucGlmer}}. There are many methods for
##' \code{\link{strucGlmer-class}} objects including
##' \code{\link{print.strucGlmer}}; \code{\link{summary}};
##' \code{\link{pars.strucGlmer}}; \code{\link{fixef.strucGlmer}};
##' \code{\link{covar.strucGlmer}}; \code{\link{loads.strucGlmer}};
##' \code{\link{ranef.strucGlmer}}; \code{\link{nobs.strucGlmer}};
##' \code{\link{residuals.strucGlmer}};
##' \code{\link{fitted.strucGlmer}}; \code{\link{vcov.strucGlmer}};
##' \code{\link{VarCorr.strucGlmer}}; \code{\link{isREML.strucGlmer}};
##' \code{\link{df.residual.strucGlmer}};
##' \code{\link{deviance.strucGlmer}};
##' \code{\link{logLik.strucGlmer}}; \code{\link{formula.strucGlmer}};
##' \code{\link{family.strucGlmer}}; \code{\link{weights.strucGlmer}};
##' \code{\link{model.matrix.strucGlmer}};
##' \code{\link{simulate.strucGlmer}}. Also see:
##' \code{\link{getStrucGlmerPar}}; \code{\link{parPerTerm}};
##' \code{\link{covarPerTerm}}; \code{\link{loadsPerTerm}};
##' \code{\link{nRePerTrm}}; \code{\link{reIndsPerTrm}};
##' \code{\link{getOffset}}; \code{\link{simStrucParsedForm}};
##' \code{\link{printReTrm}}; \code{\link{getReTrm}};
##' \code{\link{compressStrucGlmer}}.}
##'
##' \item{random-effects structures}{One may specify new special
##' random effects structures by defining new methods for the
##' \code{\link{setReTrm}} generic. To get started constructing such
##' methods see \code{showSkeleton("setReTrm")}. The
##' \code{\link{findReTrmClasses}} function may be used to list all
##' methods that are currently available. Also see:
##' \code{\link{mkReTrmStructs}}; \code{\link{getModMatAndGrpFac}};
##' \code{\link{packReTrm}}; \code{\link{update.reTrmStruct}};
##' \code{\link{setLowerDefault}}; \code{\link{simReTrm}};
##' \code{\link{mkSparseTrans}}.}
##'
##' \item{Structured sparse matrices}{Many contemporary mixed models
##' and latent variable models involve large matrices, which at each
##' iteration of an optimization algorithm need to (1) have their
##' values updated and (2) be involved in matrix algebra. Fortunately
##' these large matrices have structure that can be exploited for
##' developing better fitting algorithms. The first type of structure
##' is that these matrices are constructed by applying a small number
##' of simple operations on relatively small building-block
##' matrices. This basis in build-block matrices makes it simpler to
##' update mixed model matrices. The second type of structure is that
##' these mixed model matrices are often sparse, in that they contain
##' a relatively large number of zeros. This sparsity allows for more
##' efficient matrix operations involved in the fitting
##' proceedure. There is quite alot of machinery in the \code{Matrix}
##' package for exploiting sparsity. There is less machinery for
##' exploiting the building-block structure of mixed-model
##' matrices. The \code{\link{strucMatrix}} objects in the
##' \code{lme4ord} package provide such machinery.}
##'
##' \item{utilities}{A utility is a function that does _not_ make use
##' of functions from any other module. The purpose of separating
##' utilities in this way is to preserve modularity. Currently, the list of
##' utilities is: \code{\link{mkParInds}}; \code{\link{getInit}};
##' \code{\link{getInit}}; \code{\link{factors}},
##' \code{\link{stanCov}}; \code{\link{nChoose2Inv}};
##' \code{\link{familySimFun}}; \code{\link{triInds}};
##' \code{\link{countUnique}}; \code{\link{flattenIntVec}};
##' \code{\link{assignWith}}; \code{\link{showSkeleton}};
##' \code{\link{listTranspose}}; \code{\link{denseDiagInds}};
##' \code{\link{edgeTipIndicator}}; \code{\link{covExpDecay}};
##' \code{\link{orthProcrustesRotMat}}.}
##'
##' }
##' 
##' @docType package
##' @name lme4ord
##' @aliases lme4ord package-lme4ord lme4ord-package
NULL

