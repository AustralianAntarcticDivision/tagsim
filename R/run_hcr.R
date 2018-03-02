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
  storage <- create_storage_hcr(control, n_reps)
  ## loop over the simulations
  for(i in 1:n_reps){
    #print(paste("## Simulation ",iSim," ##",sep=""))
    ##* the first index is already created
    model <- create_model_hcr(control)
    ### Annual loop over the years
    for(y in 2:(control[["n_years"]]+1)) {
      #print(paste("## Annual loop:  Year ",control$years[y]," ##",sep=""))
      ## Simulation Order
      ### 1 New season processes
      ## 1.1 create a temp vector to store the populationnumbers
      N_yr <- model$N[y-1]
      ## 1.2 Estimate recruitment from last season numbers (there is no growth)
      rec_yr <- est_recruits(type=control[["rec_pars"]]$type, N_area = N_yr,
                             rec_pars=control[["rec_pars"]])
      ## there is no movement and only a single area so add the recruits
      ## 1.3 add the recrutiment
      ##* consider ssaving the recruitments
      N_yr <- N_yr + rec_yr
      ## 2 Harvest and natural mortality
      ## 2.2 calculate & remove total catch & apply natural mortality
      ## apply 1/2 M, the F, then 1/2 M
      N_yr <- N_yr * exp(-0.5*control[["pop_pars"]]$nat_mort)
      ## calculate then remove the catch
      catch_yr <- calc_catch(control, N_yr, model$abund_est[y-1])
      N_yr <- N_yr - catch_yr
      ## apply 1/2 natural mortality
      final_N <- N_yr * exp(-0.5*control[["pop_pars"]]$nat_mort)
      ## 3 Assign the end of season population size to the storage
      model$N[y] <- final_N
      model$catch[y] <- catch_yr
      model$recruits[y] <- rec_yr
      ## note the assessment happens at the end of the year (after both F&M)
      model$abund_est[y] <- do_assessment(control, model, year=y)
    }
    storage <- store_rep_hcr(storage, control, model, rep=i)
  }
  ## return the storage
  storage
}
