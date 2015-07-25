##' @aliases <%= cls %>
##' @inheritParams setReTrm
##' @family setReTrm
##' @section Mixed model formula usage:
##' <%= form %>
##' \describe{
##' <%= paste("\\item{", eval(parse(text = arg)), ":}{", eval(parse(text = desc)), "}", sep = "") %>
##' }
##' @section Parameters:
##' \describe{
##' <%= paste("\\item{Covariance:}{", eval(parse(text = covarDesc)), "}", sep = "") %>
##' <%= paste("\\item{Loadings:}{", eval(parse(text = loadsDesc)), "}", sep = "") %>
##' }
