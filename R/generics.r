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

#' Calculate recrutiment from previous years population
#'
#' Function returns number of recruits given last year's population size
#' carrying capacity K, natural survivorship S and resilience r.
#' rk is the number of recruits per unit spawner when B = K.
#' var controls whether the recruitment is stochastic
#' (var = F is used when calculating recruitment parameters)
#' @param type recruitment function either "constant", "logistic", "bevholt" (Beverton Holt)
#' or "ricker". The "lognorm" method is independent of population size
#' @param N_area population size by region
#' @param rec_pars list of recruitment parameters
#' @param var is recruitment stochastic (default=FALSE)
#' @export
est_recruits = function(type, N_area, rec_pars, var=FALSE){
  ## add some checks
  N <- sum(N_area)
  ## calculate recruitment based on type
  switch(type,
         constant = recruits <- rec_pars[["initial"]],
         logistic = {
           recruits <- N*rec_pars[["rk"]]*(1 + rec_pars[["resilience"]]*(1 - N/rec_pars[["K"]]))
         },
         bevholt = {
           bta <- rec_pars[["resilience"]]*rec_pars[["rk"]]
           b <- 1./(bta/(rec_pars[["rk"]]*rec_pars[["K"]]) - 1./rec_pars[["K"]])
           a <- b*bta
           recruits<-a*N/(b + N)
         },
         ricker = {
           a <- rec_pars[["resilience"]]*rec_pars[["rk"]]
           b <- log(a*rec_pars[["K"]]/(rec_pars[["rk"]]*rec_pars[["K"]]))/rec_pars[["K"]]
           recruits<-a*N*exp(-b*N)
         }
  )
  ## prevent negative recruitment
  recruits <- ifelse(recruits>0, recruits, 0)
  ## Include variability when required
  if(var) {
    recruits <- recruits * rlnorm(1,rec_pars[["mu"]],rec_pars[["s"]])
  }
  ## return the recruits
  recruits
}











