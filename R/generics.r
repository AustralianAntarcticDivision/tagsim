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












