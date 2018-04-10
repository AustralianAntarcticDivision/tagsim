## functions to implement harvest strategies

#' Select fishing locations
#'
#' Fishing locations based on a strategy and proportion area fished
#' @param type fishing strategy
#' @param prop proportion of cells to be fished
#' @param regions vector of regions
#' @param revisit_prop proportion of locations to revisit
#' @param old_locs vector of locations fished previous year
#' @export
fishing_locs <- function(type, regions, prop, revisit_prop, old_locs){
  ## checks
  if(floor(sum(prop)*length(regions)) <= 0) stop("must choose prop so at least one region is selected")
  ## define the number of hauls, rounding down
  n_hauls <- floor(prop * length(regions))
  ## random fishing
  if(type=="random"){
    fished <- sort(sample(x=regions, size=n_hauls, replace=FALSE))
    ## create a vector of fishing locations
    locs <- rep(0, length(regions))
    ## set the fished locations to 1
    locs[fished] <- 1
    ## normalise
    locs <- locs/sum(locs)
  }  else if(type=="user"){
    locs <- prop
  }else     stop("not yet implemented")
  #else if(type=="revisit"){
    ## checks
#     if(any(!old_locs %in% regions)) stop("old_locs must be a subset of regions")
#     #if((prop + revisit_prop) > 1) stop("prop + revisit_prop must be <= 1")
#     ## calc the revisit proportions
#     revisit_hauls <- floor(revisit_prop*n_hauls)
#     revisit_locs <- sample(old_locs, size=revisit_hauls, replace=FALSE)
#     new_hauls <- n_hauls-revisit_hauls
#     unfished_locs <- which(!regions %in% old_locs)
#     new_locs <- sample(unfished_locs, new_hauls, replace=FALSE)
#     ## there are going to be instances where there are less new_locs
#     fished <- sort(c(revisit_locs, new_locs))
  #}else if(type=="target"){## this strategy revisits the high catch cells
    ##this also needs two additional arguments
    ## perhaps the best way to implement this is directly specify
    ## which cells to revisit
  ## return vector of locations with fished locs = 1
  locs
}

#' Calculate catch
#'
#' Calculate catch numbers for a year
#' @param control a control file
#' @param N_area true population numbers, can be by area
#' @param est_area estimated population numbers, can be by area
#' @export
calc_catch <- function(control, N_area, est_area){
  ## total population numbers and estimate
  N <- sum(N_area)
  est <- sum(est_area)
  ## define the maximum catch
  max_catch <- N * control[["harvest_pars"]]$max
  ## define catch
  catch <- NULL
  ## calc catch based on type
  switch(control[["harvest_pars"]]$type,
         const_exploit = catch <- est * control[["harvest_pars"]]$exploit_rate,
         TAC = catch <- control[["harvest_pars"]]$exploit_rate
  )
  ## restrict catch to max harvest rate
  catch <- min(max_catch, catch)
  ## return the catch
  catch
}

#' Calculate catch by area
#'
#' Calculate the catches for each spatial cell
#' @param method either "Uniform" or in proporiton to population "Density"
#' @param catch total catch to be taken
#' @param pop vector of population numbers by area
#' @param locs vector of locations to extract catch from (zeroes denote no catch taken from area)
#' @param max_exploit maximum harvest rate (can be a vector for certain methods)
#' @param harvest if some cells will be depleted below max_exploit will the residual catch be taken elsewhere
# calc_catch <- function(method, catch, pop, locs, max_exploit, harvest=FALSE){
#   ## think about additional checks
#   ##if(abs(sum(locs) - 1) < 0.0001) stop("locs must sum to ~ 1")
#   ##if(catch > sum(pop*max_exploit)) stop("catch is greater than available population for value of max_exploit")
#   ## define the max catch by area
#   max_catch <- round(pop * max_exploit)
#   ## calculate the raw catches
#   if(method=="uniform"){
#     raw_catch <- round(catch * (locs / sum(locs)))
#   }else if(method=="density"){
#     ## for density
#     density_fished <- locs * pop / sum(locs*pop)
#     ## raw catch is in proportion to density
#     raw_catch <- round(density_fished * catch)
#   }
#   ## check the catch is available
#   check_catch <- max_catch - raw_catch
#   if(any(check_catch < 0)){
#     ## if catch is harvested then assign overcatch and recheck
#     if(harvest){
#       stop("reallocation of overcatch not yet defined")
#       ## define the overcatch
#       #       overcatch <- -1 * sum(check_catch[which(check_catch<0)])
#       #       while(overcatch>0){
#       #         ## identify locs with spare catch
#       #         spare_catch_locs <- which(check_catch > 0)
#       #         ## rescale the proportions, set other props to zero and
#       #         ##spare_locs <- locs[!spare_catch_locs] <- 0
#       #       }
#       #       catch <- obj
#     }else if(!harvest){
#       ## harvest==FALSE -> excess catch is not harvested
#       catch <- pmin(raw_catch, max_catch)
#     }
#   }else{
#     ## no overcatch raw catch is harvested
#     catch <- raw_catch
#   }
#   ## return a vector of catch numbers by area
#   return(catch)
# }
#

