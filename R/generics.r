## Functions written by Bill de la Mare from Rpopsim

#' Baranov catch equation
#'
#' Returns catch for given F for population size N
#' with natural and fishing mortality M and F respectively
#' @param N population size
#' @param M natural mortality (instantaneous)
#' @param F fishing mortality (instantaneous)
#' @export
baranov<-function(N,M,F)
{
  N*F/(M+F)*(1-exp(-M-F))
}

## ************************************************ ##
## New S3 Methods

#' Bootstrap method for single release tag-return studies
#'
#' Bootstrap method for single release tag-return studies
#'
#' Estimate confidence intervals using a non-parametric bootstrap. This method
#' incoporates uncertainty in tag-induced mortality, natural mortality and tag
#' shedding as a series of Bernoulli trials. @seealso \code{\link{single_release}}
#' for more details.
#' @param x an object of class srelease
#' @param nboot number of bootstrap samples
#' @param ... additional parameters
#' @aliases bootstrap.srelease
#' @export
bootstrap <- function(x, nboot, ...)
  UseMethod("bootstrap")
