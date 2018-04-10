#' Survey a population
#'
#' Wrapper for functions to obtain estimates of population size
#' during a simulation
#' @param control a control file
#' @param model a tempory file created by \code{\link{run_sim}}
#' @param year survey year
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
             est <- single_release(tags = sum(model$releases[year-1,]),
                                   catch = round(sum(model$catch[year,])),
                                   recaps = sum(model$recaps[year,]),
                                   method = control[["assess_pars"]]$method,
                                   unit = control[["assess_pars"]]$unit,
                                   type = control[["assess_pars"]]$ricker,
                                   tag_mort = control[["assess_pars"]]$tag_mort[year],
                                   reporting = control[["assess_pars"]]$reporting[year],
                                   nat_mort = control[["assess_pars"]]$nat_mort[year],
                                   chronic_shed = control[["assess_pars"]]$chronic_shed[year])
             
           }
         },
         multi_tag = {
           if(year==1){
             est <- 0
           }else if (year==2){
             est <- single_release(tags = sum(model$releases[year-1,]),
                                   catch = round(sum(model$catch[year,])),
                                   recaps = sum(model$recaps[,year]),
                                   method = control[["assess_pars"]]$method,
                                   unit = control[["assess_pars"]]$unit,
                                   type = control[["assess_pars"]]$ricker,
                                   tag_mort = control[["assess_pars"]]$tag_mort[year],
                                   reporting = control[["assess_pars"]]$reporting[year],
                                   nat_mort = control[["assess_pars"]]$nat_mort[year],
                                   chronic_shed = control[["assess_pars"]]$chronic_shed[year])

           }else if (year>2){
             # must reduce the size of the release and recapture matrix to exclude 
             # years for which no data have been simulated yet
             # with releases we only want those from years prior to the current year 
             # as we dont count within season recaps so we will only have recaps from 1:current year-1
             # catch from the second year onwards is only relevant too so we remove catch from year 1 
             est <- multi_release(tags = cbind(model$releases,model$recaps)[1:(year-1),1:(year+1)],
                                  hauls = model$catch[2:year],
                                  pars=control[["assess_pars"]])
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
