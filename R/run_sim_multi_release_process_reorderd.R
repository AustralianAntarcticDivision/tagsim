#' Run a simulation model
#'
#' @param control control file
#' @param n_reps number of MCMC replicates
#' @param run_assessment TRUE if you want the simulation to estimate biomass and FALSE if estimates are not required
#' @importFrom stats rbinom
#' @export
#' 
run_sim_multi_release_new_process_order <- function(control, n_reps, run_assessment){
  ## storage for the final results (do we want this by region?)
  storage <- create_storage(control, n_reps)
  ## loop over the simulations
  for(i in 1:n_reps){
    ## the first index is already created
    model <- create_model_multi_release(control)
    ### Annual loop over the years
    for(y in 1:(control[["n_years"]])) {
      ## Simulation Order
      
      ## 1.1 create a temp vector to store the population and tag numbers in year one this is the initial pop size in the control file
      if(y==1){
        # N_true vector starts with the initial biomass
        temp_pop <- model$N_true[y,]
      }else{
        
        temp_pop <- model$N_end_season[y-1,]
       
      }
      
      # multi_release function will be calculating tags available for each release event
      # tags available is a matrix in control file that in n x n where n = years of the survey 
      avail_tags <- model$tags_available
      
      ## 1.2 Move untagged  and tagged population 
      # temp_untagged   <- ceiling(move_N(temp_untagged, control[["movement"]]))
      temp_pop   <- move_N(temp_pop, control[["movement"]])
      
      
      # 1.2 Estimate recruitment
      if(y==1){
        rec <-est_recruits(type=control[["rec_pars"]]$type,
                           rec_pars=control[["rec_pars"]],
                           var=control[["rec_pars"]]$variation)}
      else{
        # replace individuals taken out as recruits
        # rec <- model$N_true[y-1]*log(exp(control[["pop_pars"]]$nat_mort)) 
        rec <- model$N_true[y-1]-temp_pop
      }
      
      
      ## assign it to areas (this can be replaced with a function)
      # rec_area <- ceiling(rec * control$rec_area)
      rec_area <- rec * control$rec_area
      
      # 1.3 add the recrutiment (recruits don't move)
      temp_pop <- temp_pop + rec_area
      
      ## store true total population to be compared with assessment estimates  
      true_N <- temp_pop
      
      ### New season processes
      ## 2 Harvest 
      ## 2.1 extract last seasons abundance estimate
      if(y>2){
        switch(control[["harvest_pars"]]$pop,
               est = est_abund <- model$abund_est[y-1,],
               true = est_abund <- true_N) 
      }else{
        est_abund <- true_N
      }
      ## 2.2 caculate catch based on harvest rate/strategy 
      temp_catch <- calc_catch(control, temp_pop, est_abund)
      ## 2.3 assign catch removal by area
      catch_by_area <- temp_catch/control[["regions"]]
      
      ## 2.4 release tags make this stochastic 
      if(control[["tag_pars"]]$catch_dependent==TRUE){
        tag_releases <- round(catch_by_area * control[["harvest_pars"]]$tag_rate)
      }else{
        tag_releases <- round(true_N*control[["harvest_pars"]]$tag_rate)
      }
      # recaps by year instead of area
      # recaps_by_yr <- rep(0,y)
      # 2.5 recapture tags based on tagged fish available from the previous season
      # 2.5.1 create storage for recaptures from each cohort of releases 
      recaps_by_yr <- rep(0, control[["n_years"]])
      
      # 2.5.2 calculate the probability of recapture based on available tags from each cohort in
      # the previous season (for year 1 this is zero) and the expected number of recaptures 
      if(y==1){
        rprob_by_yr <- avail_tags[,y] / true_N
        expected_recaps <- 0
      }else{
        rprob_by_yr <- avail_tags[,y-1] / true_N
        expected_recaps <- sum(avail_tags[,y-1])*control[["harvest_pars"]]$exploit_rate
      }
      
      for(k in 1:length(recaps_by_yr)){
        recaps_by_yr[k] <- rbinom(n=1, size=round(catch_by_area),
                                  prob=rprob_by_yr[k])
      }
      
      if(any(is.infinite(rprob_by_yr))) stop("Inf probability of recapture")
      # 2.6 in the first year remove catch from the untagged population and add releasese to the tagged population 
      # and add tagged and untagged to get a total end of season population 
      # from year 2 onwards recaptures are added back into the untagged population and removed from the tagged 
      # population
      
      if(y==1){
        # no recaptures in year 1
        temp_untagged <- temp_pop - catch_by_area
        
        # temp_tagged <- tag_releases
        
      }else{
        temp_untagged <- temp_pop - catch_by_area + sum(recaps_by_yr)
        # account for the tagged population - this is done in the tags available below
        # temp_tagged <- tag_releases - sum(recaps_by_yr)
      }
      
      # account for the untagged population
      # # add releases back in to the total population
      # total_pop <- temp_untagged + temp_tagged
      # 
      # 
      # 1.4 apply all natural mortality - check that you are not applying this twice in the first year 
      temp_untagged<- temp_untagged* exp(-control[["pop_pars"]]$nat_mort)
      
      # 3. Tags available 
      # caluclate the tags that will be available for recapture in the following season 
      ## calculate the available tags based on Ricker fishery type
      ## calculate the available tags for each cohort (row)
      # for year 1 just use tag releases and for for all years beyone this then 
      # use the tags available from the previous season from each cohort
      if(y==1){
        ## tag-induced mortality is applied in the release year only
        temp_tags <- tag_releases * exp(-control[["tag_pars"]]$tag_mort[y])
        ## remove recaptures scaled for reporting rate (could improve)
        temp_tags <- temp_tags - (recaps_by_yr[y] / control[["tag_pars"]]$reporting[y])
        ## apply natural mortality etc
        avail_tags[y,1] <- temp_tags*exp(-control[["tag_pars"]]$nat_mort[y]
                                         -control[["tag_pars"]]$chronic_shed[y]
                                         -control[["tag_pars"]]$chronic_mort[y])
      }else{
        temp_tags <- avail_tags[,y-1]
        
        for (j in 1:length(temp_tags)){
          
          if(j<y){
            ## repeat the removals above
            temp_tags_discount <- temp_tags[j] - (recaps_by_yr[j] / control[["tag_pars"]]$reporting[y])
            ## apply natural mortality etc
            avail_tags[j,y] <- temp_tags_discount*exp(-control[["tag_pars"]]$nat_mort[y]
                                                      -control[["tag_pars"]]$chronic_shed[y]
                                                      -control[["tag_pars"]]$chronic_mort[y])
          }else if(j==y){
            ## tag-induced mortality is applied in the release year only
            temp_tags_discount <- tag_releases * exp(-control[["tag_pars"]]$tag_mort[y])
            ## remove recaptures scaled for reporting rate (could improve)
            temp_tags_discount <- temp_tags_discount - (recaps_by_yr[y] / control[["tag_pars"]]$reporting[y])
            ## apply natural mortality etc
            avail_tags[j,y] <- temp_tags_discount*exp(-control[["tag_pars"]]$nat_mort[y]
                                                      -control[["tag_pars"]]$chronic_shed[y]
                                                      -control[["tag_pars"]]$chronic_mort[y])
          }else if(j>y){
            
            # do nothing   
            
          }
        }
      }
      
      # add tagged population 
      total_pop <- temp_untagged + sum(avail_tags[,y])
      
      
      # ## 4 Assign the end of season population size to the storage
      model$N_end_season[y,] <- total_pop
      model$N_true[y,] <- true_N
      model$releases[y,] <- tag_releases
      # model$ recaps is a 4 dim array with year of release as rows year of recapture as columns and the release region and recapture region in the 3rd and fourth dim of the array
      # if(y==1){model$tags_available[y,] <- temp_tags}
      # this should update the matrix each year
      model$tags_available <- avail_tags
      model$tags_available_store[y,]<- sum(avail_tags[,y])
      model$recaps[,y] <- recaps_by_yr
      model$recaps_store[y,] <- sum(recaps_by_yr)
      if(control[["assess_pars"]]$method=="Chapman" & control[["assess_pars"]]$unit %in% c("kg","tonnes")){
        model$catch[y,] <- catch_by_area*control[["assess_pars"]]$mean_wt
      }else{model$catch[y,] <- catch_by_area}
      model$recruits[y,] <- rec
      model$expected_recaps[y,] <- round(expected_recaps,0)
      ## do the assessment
      ## 5 Estimate abundance and store the result
      if (run_assessment==TRUE & y>=2){
        temp_abund <- do_assessment(control, model, year=y)
        # tags_available <- temp_abund$Avail_tags
        if (y==2){
          temp_abund <- temp_abund$Estimate[["N_hat"]]
        }else{
          temp_abund <- temp_abund$Est[["Combined"]]
        }
        ## implement this properly
        if(sum(temp_abund) == 0|is.na(temp_abund)){
          final_abund <- 0
        }else{
          final_abund <- temp_abund
        }
        ## add the abundance estimates
        model$abund_est[y,] <- final_abund
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
