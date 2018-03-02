#' Create a simplified control file for a single-area population without tags
#'
#' The control file specifies the simulation parameters
#'
#' A simulation is run using the function (see \code{\link{run_hcr}})
#' @param years vector of years (integers starting with 1)
#' @param pop_pars list of population parameters containing the initial population
#' size ("initial") and instantaneous natural mortality ("nat_mort")
#' @param rec_pars list of recruitment parameters containing the name of the mechanism
#' for calculating recruitment ("type"; can be either "constant", "logistic", "bevholt"
#' or "ricker") see \code{\link{est_recruits}}
#' @param harvest_pars harvest parameters a list containing "type" ("const_exploit" =
#' constant exploitation rate, "TAC" = constant catch, ...),  (see \code{\link{calc_catch}})
#' @param assess_pars list of assessment parameters (see \code{\link{do_assessment}})
#' @param stoch_rec stochasitc recruitment (default=TRUE)
#' @aliases control
#' @export
create_control_hcr <- function(years,
                               pop_pars,
                               rec_pars,
                               harvest_pars,
                               assess_pars,
                               stoch_rec = TRUE){
  ## add additional error checking
  ## create list object
  control <- list("years" = years,
                  "pop_pars" = pop_pars,
                  "rec_pars" = rec_pars,
                  "harvest_pars" =harvest_pars,
                  "assess_pars" = assess_pars,
                  "n_years" = length(years))
  ## add a class
  class(control) <- "control"
  ## return the control file
  control
}

#' Create a simplified model object for a single-area population without tags
#'
#' Create a list of arrays of the appropriate dimensions for storing population
#' and fishery information
#'
#' Note tags could be a multi-dimensional array to account for tag movement
#' @param control a control file
#' @export
create_model_hcr <- function(control){
  ## create storage for population size, recruitment and assessment
  init_N <- rep(0, control[["n_years"]] + 1)
  init_R <- rep(0, control[["n_years"]] + 1)
  init_A <- rep(0, control[["n_years"]] +1)
  ## fill year zero
  if(control[["rec_pars"]]$stochastic_rec){
    init_N[1] <- control[["pop_pars"]]$initial *
                  rlnorm(1, control[["rec_pars"]]$mu, control[["rec_pars"]]$s)
  }else{
    init_N[1] <- control[["pop_pars"]]$initial
  }
  ## add initial recruitment (not getting used currently)
  init_R[1] <- est_recruits(type=control[["rec_pars"]]$type,
                            N_area = 1, rec_pars=control[["rec_pars"]])
  ## initial assessment knows pop size without error ** can change this
  init_A[1] <- init_N[1]
  ## create the object
  obj <- list("N" = init_N,
              "recruits" = init_R,
              "catch" = rep(0, control[["n_years"]] + 1),
              "abund_est" = init_A
  )
  ## return the object
  return(obj)
}

#' Create an object to store a hcr simulation results
#'
#' Create an object to storea hcr simulation results
#' @param control a control file
#' @param n_reps number of replicate simulations
#' @export
create_storage_hcr <- function(control, n_reps){
  ## create a list to store the simulation output
  obj<- list("true_N" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
             "est_N" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
             "recruits" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
             "catch" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps))
  ## return the object
  obj
}

## not going to be very memory efficient
store_rep_hcr <- function(storage, control, model, rep){
  ## store true and estimated N, recruits and catch
  storage$true_N[,rep] <- model[["N"]]
  storage$est_N[,rep] <- model[["abund_est"]]
  storage$recruits[,rep] <- model[["recruits"]]
  storage$catch[,rep] <- model[["catch"]]
  ## return the object
  storage
}

