## tag-based estimators of population size 

#' Petersen estimators
#' 
#' Single tag recapture estimates Petersen, Lincoln Petersen estimators
#' @param tags number of marked animals released (note these must be whole 
#' numbers for the bootstrap to work)
#' @param catch vector of the number or weight of animals captured and checked 
#' for tags by haul on the second survey
#' @param recaps vector of number of marked animals recaptured by haul
#' @param mean_wt mean weight of a single fish (optional argument for Chapman
#'  weight method)
#' @param check_type of input checks "srelease" (default= single release checks)
#' other options are "mrelease" (multiple release checks) and "zero_recaps"
#' (include 
#' @name petersen
NULL
#> NULL


#' @export
#' @rdname petersen
chapman_n <- function(tags, catch, recaps, check_type="srelease"){
  ## check function inputs
  switch(check_type,
         srelease = {check <- check_srelease_inputs(tags, catch, recaps)},
         mrelease = {check <- check_mrelease_inputs(tags, catch, recaps)})
  ## if check passes calculate population size
  if(check){
    ## calculate the estimate
    N_hat <- (((tags + 1)*(catch + 1))/(recaps + 1)) - 1 
    ## calculate the variance
    var_N <- ((tags + 1)*(catch + 1)*(tags - recaps)*(catch - recaps))/
      (((recaps + 1)^2)*(recaps + 2))
    obj <- c(N_hat, var_N)
  }else{
    obj <- c(NA, NA)
  }
  names(obj) <- c("N_hat", "var_N")
  obj
}

#' @export
#' @rdname petersen
chapman_wt <- function(tags, catch, recaps, mean_wt=0, check_type="srelease"){
  ## check function inputs
  switch(check_type,
         srelease = {check <- check_srelease_inputs(tags, catch, recaps)},
         mrelease = {check <- check_mrelease_inputs(tags, catch, recaps)})
  ## if check passes calculate population size
  if(check){
    ## calculate the estimate
    N_hat <- (((tags + 1)*(catch + mean_wt))/(recaps + 1)) - mean_wt
    ## calculate the variance
    var_N <- ((tags + 1)*(catch + mean_wt)*(tags - recaps)*(catch - mean_wt*recaps))/
      (((recaps + 1)^2)*(recaps + 2))
    obj <- c(N_hat, var_N)
  }else{
    obj <- c(NA, NA)
  }
  names(obj) <- c("N_hat", "var_N")
  obj
}

#' @export
#' @rdname petersen
petersen <- function(tags, catch, recaps, check_type="srelease"){
  ## check function inputs
  switch(check_type,
         srelease = {check <- check_srelease_inputs(tags, catch, recaps)},
         mrelease = {check <- check_mrelease_inputs(tags, catch, recaps)})
  ## if check passes calculate population size
  if(check){
    ## calculate the estimate
    N_hat <- tags * catch / recaps 
    ## variance from Ricker 1975
    var_N <- (tags^2) * catch * (catch - recaps) / recaps^3 
    ## return the estimate and variance
    obj <- c(N_hat, var_N)
  }else{
    obj <- c(NA, NA)
  }
  names(obj) <- c("N_hat", "var_N")
  obj
}
