## simplified simulation loop to demonstrate harvest control rules

#' Run a single-area simulation with a harvest control rule
#'
#' Run a simulation using a harvest control rule to modify catch during the
#' simulation
#' @param control control file
#' @param n_reps number of replicate simulations
#' @param save save the output (default=FALSE)
#' @param path path for the output
#' @param filename filename for the control file
#' @export
#' @examples
#' ## add an example
run_hcr <- function(control, n_reps, save=FALSE, path=NULL, filename=NULL){
  ##* storage for the final results
  storage <- create_storage(control, n_reps)
  ## loop over the simulations
  for(i in 1:n_reps){
    #print(paste("## Simulation ",iSim," ##",sep=""))
    ##* the first index is already created
    model <- create_model(control)
    ### Annual loop over the years
    for(y in 2:(control[["n_years"]]+1)) {
      #print(paste("## Annual loop:  Year ",control$years[y]," ##",sep=""))
      ## Simulation Order
      ### 1 New season processes
      ## 1.1 create a temp vector to store the populationnumbers
      temp_N <- model$N[y-1,]
      ## 1.2 Estimate recruitment from last season numbers (there is no growth)
      rec <- est_recruits(type=control[["rec_pars"]]$type,
                          N_area = temp_N,
                          rec_pars=control[["rec_pars"]])
      ## there is no movement and only a single area so add the recruits
      ## 1.3 add the recrutiment
      ##* consider ssaving the recruitments
      temp_N <- temp_N + rec
      ## 2 Harvest and natural mortality
      ## 2.1.1 extract last seasons abundance estimate
      est_abund <- model$abund_est[y-1,]
      ## 2.1.2 update the harvest rate/strategy based
      ##* check that the maximum catch is implemented here
      catch_limit <- calc_catch(control, temp_N, est_abund)
      ## 2.2 calculate & remove total catch & apply natural mortality
      if(control[["harvest_pars"]]$ricker ==1){
        ## ricker type 1 all fishing first then natural mortality
        temp_N <- temp_N - catch_limit
        ## apply natural mortality
        temp_N <- temp_N * exp(-control[["pop_pars"]]$nat_mort)
        ##cat("final_N ", final_N, "\n")
        }else if(control[["harvest_pars"]]$ricker ==2){
          ## ricker type 2, 1/2 M, the F, then 1/2 M
          ## apply 1/2 natural mortality
          temp_N <- temp_N * exp(-0.5*control[["pop_pars"]]$nat_mort)
          ## remove catch
          temp_N <- temp_N - catch_limit
          ## apply 1/2 natural mortality
          temp_N <- temp_N * exp(-0.5*control[["pop_pars"]]$nat_mort)
          ##cat("final_N ", final_N, "\n")
      }
      ## 3 Assign the end of season population size to the storage
      model$N[y,] <- final_N
      model$catch[y,] <- catch_by_area
      model$recruits[y,] <- rec_area
      ## do the assessment
      ##* modify the assessment so it can't go below zero
      model$abund_est[y,] <- do_assessment(control, model, year=y)
    }
    storage <- store_rep(storage, control, model, rep=i)
  }
  ## return the storage
  storage
}
