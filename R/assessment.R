#' Survey a population
#'
#' Wrapper for functions to obtain estimates of population size
#' during a simulation
#' @param control a control file
#' @param model a tempory file created by \code{\link{run_sim}}
#' @param year survey year
#' @import tagr
#' @importFrom stats rnorm 
#' @export
do_assessment <- function(control, model, year){
  ## some checks
  ## could calculate the true (stratified) population size
  switch(control[["assess_pars"]]$type,
         single_tag = {
           ## Petersen estimate of abundance
           if(year==1){
             est <- 0
           }else{
             est <- tagr::single_release(tags = sum(model$releases[year-1,]),
                                         catch = round(sum(model$catch[year,])),
                                         recaps = sum(model$recaps[year,]),
                                         method = "Petersen",
                                         unit = "numbers",
                                         type = control[["harvest_pars"]]$ricker,
                                         tag_mort = control[["tag_pars"]]$mort,
                                         reporting = control[["tag_pars"]]$report,
                                         nat_mort = control[["pop_pars"]]$nat_mort,
                                         chronic_shed = control[["tag_pars"]]$shed)

           }
         },
         multi_tag = {
           if(year==1){
             est <- 0
           }else{
             est <- tagr::multi_release(tags = sum(model$releases[year-1,]),
                                        catch = round(sum(model$catch[year,])),
                                        recaps = sum(model$recaps[year,]),
                                        method = "Petersen",
                                        unit = "numbers",
                                        type = control[["harvest_pars"]]$ricker,
                                        tag_mort = control[["tag_pars"]]$mort,
                                        reporting = control[["tag_pars"]]$report,
                                        nat_mort = control[["pop_pars"]]$nat_mort,
                                        chronic_shed = control[["tag_pars"]]$shed)
           }
         }
         ,
         strat_tag = {
           ## this is Andrew's method
           est <- NULL
         },
         survey = {
           ## estimate biomass with some error
           pop_size <- sum(model$N[year,])
           rand_norm <- max(rnorm(1, 1, control[["assess_pars"]]$cv), 0)
           est <- pop_size * rand_norm
           #   },
           #   strat_survey = {
           #     ## as above but stratified
           #     est <- NULL
           #   },
           #   rel_abund = {
           #     ## some measure of relative abundace such as CPUE
           #     est <- NULL
         }
  )
  ## if a single estimate
  # if(length(est)==1){
  #   est <- rep(est / control[["n_regions"]])
  # }
  ## return the estimate
  est
}
