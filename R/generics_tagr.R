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
