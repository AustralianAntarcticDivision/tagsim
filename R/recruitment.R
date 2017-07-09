#' Estimate recruitment
#'
#' Estimate recruitment the total annual recruitment
#' @param type either "constant", "logistic", "bevholt" or "ricker"
#' @param N_area population by area
#' @param rec_pars recruitment parameters (from control file)
#' @export
est_recruits <- function(type, N_area, rec_pars){
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
