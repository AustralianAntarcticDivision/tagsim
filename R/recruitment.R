#' Calculate recrutiment from previous years population
#'
#' Function returns number of recruits given last year's population size N,
#' carrying capacity k and resilience r. rk is the number of recruits per unit
#' spawner when B = K.
#'
#' Modification of function written by Bill de la Mare
#' @param type recruitment function either "constant", "logistic", "bevholt"
#' (Beverton Holt) or "ricker". The "lognorm" method is independent of
#' population size
#' @param N_area population size by region
#' @param rec_pars list of recruitment parameters (note: extend this help)
#' @export
est_recruits <- function(type, N_area, rec_pars){
  ## add some checks
  N <- sum(N_area)
  ## calculate recruitment based on type
  switch(type,
         constant = recruits <- rec_pars[["initial"]],
         logistic = {
           recruits <- N*rec_pars[["rk"]]*(1 + rec_pars[["r"]]*(1 - N/rec_pars[["k"]]))
         },
         bevholt = {
           ## If r, k and rk are supplied we use them to calculate a and b
           if(all(c("r", "k", "rk") %in% names(rec_pars))){
           bta <- rec_pars[["r"]]*rec_pars[["rk"]]
           b <- 1./(bta/(rec_pars[["rk"]]*rec_pars[["k"]]) - 1./rec_pars[["k"]])
           a <- b*bta
           }
           recruits<-(rec_pars[["a"]]*N)/(rec_pars[["b"]] + N)
         },
         ricker = {
           ## If r, k and rk are supplied we use them to calculate a and b
           if(all(c("r", "k", "rk") %in% names(rec_pars))){
           a <- rec_pars[["r"]]*rec_pars[["rk"]]
           b <- log(a*rec_pars[["k"]]/(rec_pars[["rk"]]*rec_pars[["k"]]))/rec_pars[["k"]]
           }
           recruits<-rec_pars[["a"]]*N*exp(-rec_pars[["b"]]*N)
         }
  )
  ## prevent negative recruitment
  recruits <- ifelse(recruits>0, recruits, 0)
  ## Include variability when required
  if(rec_pars[["stochastic_rec"]]) {
    recruits <- recruits * rlnorm(1,rec_pars[["mu"]],rec_pars[["s"]])
  }
  ## return the recruits
  recruits
}

#' Assign recruitment proportions by area
#'
#' Calculate the proportions of recruitment by spatial cell
#' @param method current options include "uniform", "lognormal" and "user" defined
#' @param n_areas number of model areas
#' @param user_props user defined proportions
#' @param ... additional parameters
#' @export
recruit_area <- function(method, n_areas, user_props=NULL, ...){
  ## add some checks
  if(method=="uniform"){
    obj <- rep(1/n_areas, times=n_areas)
  }else if(method=="lognormal"){
    rs <- rlnorm(n_areas, ...)
    obj <- rs/sum(rs)
  }else if(method=="user"){
    #if(n_areas != length(user_props)) stop("Length of user_props is not equal to n_areas")
    ## add a tolerance to this
    ##if(sum(user_props) !=1) stop("user_props must sum to 1")
    obj <- user_props
  }else stop("method not implemented")
  ## return recruitment by area
  obj
}
