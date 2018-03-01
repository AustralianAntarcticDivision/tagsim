## functions to calculate initial parameters

#' Calculate the equilibrium population size
#'
#' Calculate the equilibrium population size. Can be used as an estimate of
#' k
#' @param pop_pars population parameters including the initial population size
#' "initial" and natural mortality "nat_mort"
#' @param rec_pars recruitment parameters
init_k <- function(pop_pars, rec_pars, max_reps=1000, tol=1){
  ## some checks
  if(pop_pars[["initial"]]==0)
    stop("Population starting value must be greater than zero")
  if(pop_pars[["nat_mort"]]<=0)
    stop("Natural mortality must be greater than zero")
  ## set stochastic recruitment to FALSE
  rec_pars[["stochastic_rec"]] <- FALSE
  ## set the initial population size
  N_temp <- pop_pars[["initial"]]
  ## reps counter
  count_reps <- 0
  ## loop to the equilibrium population size
  while(count_reps <= max_reps){
    ## calculate the recruitment
    current_rec <- est_recruits(type=rec_pars[["type"]], N_area=N_temp, rec_pars)
    current_N <- (N_temp+current_rec)*exp(-pop_pars[["nat_mort"]])
    if(abs(current_N - N_temp)<tol) return(current_N)
    ## increment the counter
    count_reps <- count_reps+1
    N_temp <- current_N
  }
  print("convergence was not achieved in the specified number of replicates")
}


