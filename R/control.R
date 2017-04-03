#' Create a control file
#'
#' The control file specifies the simulation parameters
#'
#' A simulation is run using the function (see \code{\link{run_sim}})
#' @param years vector of years (integers starting with 1)
#' @param regions vector of regions (integers starting with 1)
#' @param pop_pars list of population parameters containing the initial
#' population size ("initial") and instantaneous natural mortality ("nat_mort")
#' @param rec_pars list of recruitment parameters containing the name of the
#' mechanism for calculating recruitment ("type"; can be either "constant",
#' "logistic", "bevholt" or "ricker") 
#' Currently the resilience is steepness (for TOA=0.75 and for TOP=0.7), rk and K is carrying capacity
#' and rk is the number of recruits per unit spawner when B = K
#' see \code{\link{est_recruits}}, the
#' spatial distribution of the recruitment is assigment with "spat_dist" which
#' can be "uniform", "lognormal" or "user" specified
#' (see \code{\link{recruit_area}}). Recruitment The stochasisity in recruitment
#' is determined
#' @param harvest_pars harvest parameters a list containing "type"
#' ("const_exploit" = constant exploitation rate, "TAC" = constant catch, ...),
#' (see \code{\link{calc_catch}})
#' @param tag_pars list of tag related parameters (release, mort, shed, report)
#' @param fish_pars list of parameters relating to the fishing strategy. The
#' intent is for these to control the spatial dynamics of fishing and for
#' harvest_pars to determine total catch (see \code{\link{fishing_locs}})
#' @param move_pars list of movement parameters "type", "prob" (see
#'  \code{\link{move_matrix}})
#' @param assess_pars list of assessment parameters (see
#' \code{\link{do_assessment}})
#' @param stoch_rec stochasitc recruitment (default=TRUE)
#' @aliases control
#' @export
create_control <- function(years,
                           regions,
                           pop_pars,
                           rec_pars,
                           harvest_pars,
                           tag_pars,
                           fish_pars,
                           move_pars,
                           assess_pars,
                           stoch_rec = TRUE){
  ## add additional error checking
  ## create list object
  control <- list("years" = years,
                  "regions" = regions,
                  "pop_pars" = pop_pars,
                  "rec_pars" = rec_pars,
                  "harvest_pars" =harvest_pars,
                  "tag_pars" = tag_pars,
                  "fish_pars" = fish_pars,
                  "move_pars" = move_pars,
                  "assess_pars" = assess_pars,
                  "n_years" = length(years),
                  "n_regions" = length(regions),
                  "stochastic_rec" = stoch_rec,
                  ## define the movement matrix
                  "movement" = move_matrix(method=move_pars$type,
                                           n_rows=sqrt(length(regions)),
                                           n_cols=sqrt(length(regions)),
                                           prob=move_pars$prob,
                                           matrix=move_pars$matrix),
                  ## divide the recruitment evenly among the cells
                  "rec_area" = recruit_area(method=rec_pars["spat_dist"],
                                            n_areas=length(regions),
                                            user_props=rec_pars[["user_prop"]])
                  )
  ## add a class
  class(control) <- "control"
  ## return the control file
  control
}

#' Check control file
#'
#' Check the parameters in a control file
#' @param control a control file
#' @export
check_control <- function(control){
  ## run through checks & return
  NULL
}

#' Create a model object
#'
#' Create a list of arrays of the appropriate dimensions for storing population
#' and fishery information
#'
#' Note tags could be a multi-dimensional array to account for tag movement
#' @param control a control file
#' @export
create_model <- function(control){
  ## create storage for the
  init_N <- matrix(0, control[["n_years"]] + 1, control[["n_regions"]])
  init_R <- matrix(0, control[["n_years"]] + 1, control[["n_regions"]])
  init_A <- matrix(0, control[["n_years"]] +1 , control[["n_regions"]])
  ## fill year zero
  if(control[["stochastic_rec"]]){
    init_N[1,] <- control[["rec_area"]] * control[["pop_pars"]]$initial *
                  rlnorm(1, control[["rec_pars"]]$mu, control[["rec_pars"]]$s)
  }else{
    init_N[1,] <- control[["rec_area"]] * control[["pop_pars"]]$initial
  }
  ## add initial recruitment (not getting used currently)
  init_R[1,] <- control[["rec_area"]] *
    est_recruits(type=control[["rec_pars"]]$type,
                 rec_pars=control[["rec_pars"]],
                 var=control[["stochastic_rec"]])
  ## initial assessment knows pop size without error ** can change this
  init_A[1,] <- init_N[1,]
  ## create the object
  obj <- list("N" = init_N,
              "tags" =  matrix(0,control[["n_years"]] + 1,control[["n_regions"]]),
              "recaps" = matrix(0,control[["n_years"]] + 1,control[["n_regions"]]),
              "recruits" = init_R,
              "catch" = matrix(0,control[["n_years"]] + 1,control[["n_regions"]]),
              "abund_est" = init_A
  )
  ## return the object
  return(obj)
}

#' Create an object to store simulation results
#'
#' Create an object to store simulation results
#' @param control a control file
#' @param n_reps number of MCMC replicates
#' @export
create_storage <- function(control, n_reps){
  ##
  if(control[["assess_pars"]]$type %in% c("single_tag", "survey","multi_tag")){
    obj<- list("true_N" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
               "est_N" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
               "catch" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps))
  ## return the object
  return(obj)
  }else if(control[["assess_pars"]]$type %in% c("const_TAC")){
    obj<- list("true_N" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
               "catch" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps),
               "effort" = matrix(0, nrow=control[["n_years"]] + 1, ncol=n_reps)
               )
    ## return the object
    return(obj)
  }else stop("storage method not defined")
}

## not going to be very memory efficient
store_sim <- function(storage, control, model, sim){
  ##
  if(control[["assess_pars"]]$type %in% c("single_tag", "survey","multi_tag")){
    storage$true_N[,sim] <- rowSums(model[["N"]])
    storage$est_N[,sim] <- rowSums(model[["abund_est"]])
    storage$catch[,sim] <- rowSums(model[["catch"]])
    ## return the object
    return(storage)
  }else if(control[["assess_pars"]]$type %in% c("const_TAC")){
    ## Klaas to fill in
    ## return the object
    return(storage)
  }else stop("storage method not defined")
}
