#' Run a simulation model
#'
#' @param control control file
#' @param n_reps number of replicate simulations
#' @param save save the output (default=FALSE)
#' @param path path for the output
#' @param filename filename for the control file
#' @export
#' @examples
#' # add an example
run_sim <- function(control, n_reps, save=FALSE, path=NULL, filename=NULL){
  ## storage for the final results (do we want this by region?)
  storage <- create_storage(control, n_reps)
  ## loop over the simulations
  for(i in 1:n_reps){
    #print(paste("## Simulation ",iSim," ##",sep=""))
    ## the first index is already created
    model <- create_model(control)
    ### Annual loop over the years
    for(y in 2:(control[["n_years"]]+1)) {
      #print(paste("## Annual loop:  Year ",control$years[y]," ##",sep=""))
      ## Simulation Order
      ### 1 New season processes
      ## 1.1 create a temp vector to store the population and tag numbers
      temp_N <- model$N[y-1,]
      # this is zero until tags released, however, doesn't cause problems
      temp_tags <- model$tags[y-1,]
      ##cat("pre year ", temp_N, "\n")
      #' 1.2 Estimate recruitment from last season numbers (there is no growth)
      rec <- est_recruits(type=control[["rec_pars"]]$type,
                          N_area = temp_N + temp_tags,
                          rec_pars=control[["rec_pars"]])
      ## assign it to areas (this can be replaced with a function)
      rec_area <- ceiling(rec * control$rec_area)
      ## 1.2.1 Move untagged & tagged population
      ##cat("pre move ", temp_N, "\n")
      temp_N   <- ceiling(move_N(temp_N, control[["movement"]]))
      ##cat("post move ", temp_N, "\n")
      temp_tags <- ceiling(move_N(temp_tags, control[["movement"]]))
      #' 1.3.1 add the recrutiment (recruits don't move) consider saving this
      temp_N <- temp_N + rec_area
      ##cat("post recruits ", temp_N, "\n")
      ## total population
      ## 2 Harvest and natural mortality
      ## 2.1 calculate population size by strata (tagged + untagged)
      temp_pop <- temp_N + temp_tags
      ## 2.1.1 extract last seasons abundance estimate
      est_abund <- model$abund_est[y-1,]
      ## 2.1.2 update the harvest rate/strategy based
      temp_catch <- calc_catch(control, temp_pop, est_abund)
      ## specify how catch is allocated amoung regions
      if(sum(model$catch[y-1,])==0){
        previous_locs <- rep(0, control[["n_regions"]])
      }else{
        previous_locs <- model$catch[y-1,] / sum(model$catch[y-1,])
      }
      ## proportions of catch by area
      fish_by_area <- fishing_locs(type = control[["fish_pars"]]$type,
                                   regions = control[["regions"]],
                                   prop = control[["fish_pars"]]$prop,
                                   revisit_prop = control[["fish_pars"]]$revisit,
                                   old_locs = previous_locs)
      ##cat("fish_by_area", fish_by_area, "\n")
      ## 2.1.2 assign catch by area
      catch_by_area <-  fish_by_area * temp_catch
      ## also the probability of tag recapture
      rprob_by_area <- temp_tags / temp_pop
      #if(any(prob_recap < 0 | prob_recap > 1)) stop("Probability of tag recapture in cell must be between zero and one")
      ## 2.2 calculate & remove total catch & tags & apply natural mortality
      if(control[["harvest_pars"]]$ricker ==1){
        ## ricker type 1 all fishing first then
        temp_N <- temp_N - catch_by_area
        ## apply natural mortality
        temp_N <- temp_N * exp(-control[["pop_pars"]]$nat_mort)
        ## calculate tag recaptures / attrition
        ## object for recaptures
        recaps_by_area <- rep(0, control[["n_regions"]])
        for(k in 1:control$n_regions){
          recaps_by_area[k] <- rbinom(n=1, size=round(catch_by_area[k]), prob=rprob_by_area[k])
        }
        ## then remove them restricting to max number of tags
        recaps_by_area <- pmin(recaps_by_area, temp_tags)
        ## add these back to the untagged population
        ##cat("XX ", temp_N)
        temp_N <- temp_N + recaps_by_area
        temp_tags <- temp_tags - recaps_by_area
        ###*** release tags make this stochastic
        tag_releases <- round(catch_by_area * control[["tag_pars"]]$rel_rate * exp(-control[["tag_pars"]]$mort))
        ## apply natural mortality to untagged population
        ##cat("before final_N ", temp_N, "\n")
        final_N <- temp_N * exp(-control[["pop_pars"]]$nat_mort)
        ##cat("final_N ", final_N, "\n")
        ## apply natural mortality to tagged population
        final_tags <- rep(0, control[["n_regions"]])
        for(m in 1:control[["n_regions"]]){
          final_tags[m] <- rbinom(n=1, size=temp_tags[m], prob=exp(-control[["pop_pars"]]$nat_mort))
        }
      }else if(control[["harvest_pars"]]$ricker ==2){
        stop("Ricker type 2 fishery still to be implemented")
      }
      ## 3 Assign the end of season population size to the storage
      model$N[y,] <- final_N
      model$tags[y,] <- final_tags + tag_releases
      model$recaps[y,] <- recaps_by_area
      model$catch[y,] <- catch_by_area
      model$recruits[y,] <- rec_area
      ## do the assessment
      ## 4 Estimate abundance and store the result
      temp_abund <- do_assessment(control, model, year=y)
      ## implement this properly
      if(sum(temp_abund) == 0){
        final_abund <- model$abund_est[y-1,]
      }else{
        final_abund <- temp_abund
      }
      ## add the abundance estimates
      model$abund_est[y,] <- final_abund
    }
    storage <- store_rep(storage, control, model, rep=i)
  }
  ## return the storage
  storage
}
