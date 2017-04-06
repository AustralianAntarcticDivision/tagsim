#' Run a simulation model
#'
#' @param control control file
#' @param n_reps number of MCMC replicates
#' @param run_assessment TRUE if you want the simulation to estimate biomass and FALSE if estimates are not required
#' @importFrom stats rbinom
#' @export
run_sim_single_release <- function(control, n_reps, run_assessment){
  ## storage for the final results (do we want this by region?)
  storage <- create_storage(control, n_reps)
  ## loop over the simulations
  for(i in 1:n_reps){
    ## the first index is already created
    model <- create_model(control)
    ### Annual loop over the years
    for(y in 1:(control[["n_years"]])) {
      ## Simulation Order
      ### 1 New season processes
      ## 1.1 create a temp vector to store the population and tag numbers in year one this is the initial pop size in the control file
      temp_untagged <- model$N[y,]
      if(y==1){
        # zero until tags released, however, doesn't cause problems
        ## calculate the adjusted releases based on tag-induced mortality
        temp_tags <- model$tags_available[y,]
        # pop_size <- temp_untagged
      }else{
      temp_tags <- model$tags_available[y-1,] 
        # pop_size <- temp_untagged + temp_tags
        # 1.2 Estimate recruitment 
      }
      rec <-est_recruits(type=control[["rec_pars"]]$type,
                         rec_pars=control[["rec_pars"]],
                         var=control[["stochastic_rec"]])
      ## assign it to areas (this can be replaced with a function)
      rec_area <- ceiling(rec * control$rec_area)
      
      ## 1.2.1 Move untagged & tagged population
      temp_untagged   <- ceiling(move_N(temp_untagged, control[["movement"]]))
      temp_tags <- ceiling(move_N(temp_tags, control[["movement"]]))
      # 1.3.1 add the recrutiment (recruits don't move) consider saving this
      temp_untagged <- temp_untagged + rec_area
      # apply all natural mortality 
      temp_untagged <- temp_untagged* exp(-control[["pop_pars"]]$nat_mort)
      # temp_tags <- temp_tags * exp(-control[["pop_pars"]]$nat_mort/2) - removing as using tags avail as calculated in tagr which should include NM
      # maybe try the calculation of true population from here
      ## total population
      final_N <- temp_untagged
      
      ## 2 Harvest and natural mortality
      ## 2.1 calculate population size by strata (tagged + untagged)
      # temp_pop <- temp_untagged + temp_tags
      ## 2.1.1 extract last seasons abundance estimate
      if(y>2){
        est_abund <- model$abund_est[y-1,]}else{
          est_abund <- temp_untagged
        }
      ## 2.1.2 update the harvest rate/strategy based
      temp_catch <- calc_catch(control, temp_untagged, est_abund)
      ## specify how catch is allocated amoung regions
      # if(sum(model$catch[y-1,])==0){
      #   previous_locs <- rep(0, control[["n_regions"]])
      # }else{
      #   previous_locs <- model$catch[y-1,] / sum(model$catch[y-1,])
      # }
      # ## proportions of catch by area
      # fish_by_area <- fishing_locs(type = control[["fish_pars"]]$type,
      #                              regions = control[["regions"]],
      #                              prop = control[["fish_pars"]]$prop,
      #                              revisit_prop = control[["fish_pars"]]$revisit,
      #                              old_locs = previous_locs)
      ## 2.1.2 assign catch by area
      catch_by_area <- temp_catch/control[["regions"]]
      
      ###*** release tags make this stochastic
      # in the first year you apply  
      tag_releases <- round(catch_by_area * control[["tag_pars"]]$rel_rate)
      
      recaps_by_area <- rep(0, control[["n_regions"]])
      
      if(y>=2){
        rprob_by_area <- temp_tags / temp_untagged
        
        for(k in 1:control$n_regions){
          recaps_by_area[k] <- rbinom(n=1, size=round(catch_by_area[k]),
                                      prob=rprob_by_area[k])
          recaps_by_area[is.na(recaps_by_area[k])]<-0
        }
      }
      #if(any(prob_recap < 0 | prob_recap > 1))
      #stop("Probability of tag recapture in cell must be between zero and one")
      ## 2.2 calculate & remove total catch & tags & apply natural mortality
      # if(control[["harvest_pars"]]$ricker ==1){
      ## ricker type 1 all fishing first then
      ## calculate tag recaptures / attrition
      ## object for recaptures
      
      
      # ## then remove them restricting to max number of tags 
      # recaps_by_area <- pmin(recaps_by_area, temp_tags)
      ## add these back to the untagged population
      # remove catch 
      temp_untagged <- temp_untagged - catch_by_area
      # add releases back in accounting for tagging mortaility 
      temp_N <- temp_untagged + tag_releases*exp(-control[["tag_pars"]]$mort)
      # subtract recaptures from the population 
      if(y>=2){
        temp_N <- temp_N-recaps_by_area}
      # temp_tags <- temp_tags - recaps_by_area - recaptures subtraction is being handled in tagr
      # apply natural mortality to untagged population
      #cat("before final_N ", temp_untagged, "\n")
       # final_N <- temp_N
      
      # get tag releases ready for the tags_available in the following season 
      if(y==1){
        ## tag-induced mortality is applied in the release year only
        temp_tags <- tag_releases * exp(-control[["tag_pars"]]$mort)
        ## apply natural mortality etc
        temp_tags <- temp_tags*exp(-control[["pop_pars"]]$nat_mort
                                   -control[["tag_pars"]]$shed)}
      
      # ## 3 Assign the end of season population size to the storage
      model$N[y,] <- final_N
      model$releases[y,] <- tag_releases
      # model$ recaps is a 4 dim array with year of release as rows year of recapture as columns and the release region and recapture region in the 3rd and fourth dim of the array
      if(y==1){model$tags_available[y,] <- temp_tags}
      model$recaps[y,] <- recaps_by_area
      model$catch[y,] <- catch_by_area
      model$recruits[y,] <- rec_area
      ## do the assessment
      ## 4 Estimate abundance and store the result
      if (run_assessment==TRUE & y>=2){
        temp_abund <- do_assessment(control, model, year=y)
        tags_available <- temp_abund$TagsAvailable
        temp_abund <- temp_abund$Estimate[["N_hat"]]
        ## implement this properly
        if(sum(temp_abund) == 0|is.na(temp_abund)){
          final_abund <- 0
        }else{
          final_abund <- temp_abund
        }
        ## add the abundance estimates
        model$abund_est[y,] <- final_abund
        model$tags_available[y,] <- tags_available 
        # note that the model matrix stores estimates by years = rows and areas = cols
      }else{
        model$abund_est[y,] <- 0
      }
    }
    storage <- store_sim(storage, control, model, sim=i)
    # all that is currently stored is a rowsum across areas for the true N, estimated N and catch for each year of the sim 
  }
  ## return the storage
  storage
}
