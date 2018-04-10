#' Single tag release estimate of population size
#' 
#' Estimate population size from single release tag data
#' 
#' For a Type 1 fishery the adjusted releases are calculate as follows
#' 
#' 1) Initial tag-induced mortality is applied
#' 
#' 2) Within season recaptures divided by the tag reporting rate are removed
#' 
#' 3) Natural mortality and chronic tag loss and tag-induced mortality is
#'  applied.
#' 
#' For a fishery that operates throughout the year (Ricker Type 2) fishing and 
#' natural mortality competed to deplete the tagged population. The adjusted 
#' releases are calculated as follows
#' 
#' 1) Initial tag-induced mortality is applied
#' 
#' 2) Half of the within season recaptures divided by the tag reporting rate are
#' removed 
#' 
#' 3) Half of the natural mortality and chronic tag loss and tag-induced 
#' mortality is applied
#' 
#' 4) The remaining within season recaptures divided by the tag reporting rate 
#' are removed 
#' 
#' 5) The remaining natural mortality and chronic tag loss and tag-induced 
#' mortality is applied.
#' @param tags number of marked animals released (note these must be whole 
#' numbers for the bootstrap to work)
#' @param catch vector of the number or weight of animals captured and checked 
#' for tags by haul on the second survey
#' @param recaps vector of number of marked animals recaptured by haul
#' @param mean_wt mean weight of a single fish (optional argument for Chapman
#'  weight method)
#' @param prior_recaps vector of subsequent recaptures (first year this 
#' represents within season recaptures)
#' @param method method
#' @param unit measurement unit (numbers, kg, tonnes)
#' @param type do we assume fishing occurs over a short period (Ricker type 1) 
#' or is extended over the year and competes with natural mortality and tag
#' sheddding (Ricker type 2) 
#' @param tag_mort initial tag-induced mortality
#' @param reporting vector of tag reporting rate
#' @param nat_mort vector of natural mortality (instantaneous)
#' @param chronic_shed vector of chronic (ongoing) tag shedding
#' @param chronic_mort vector of chronic (ongoing) tag-induced mortality
#' @export
#' @examples 
#' ## Chapman estimates of population numbers
#' single_release(tags=100, catch=200, recaps=5, method="Chapman",
#'                type=1, unit="numbers")
#' single_release(tags=100, catch=200, recaps=5, method="Chapman", 
#'                type=1, unit="numbers", tag_mort=0.1)
#' single_release(tags=100, catch=200, recaps=5, method="Chapman", type=1,
#'                unit="numbers", tag_mort=0.1, reporting=0.8) 
#' 
#' ## Chapman estimate of the population size of carp in a lake
#' ## Load the carphauls data
#' attach(carphauls)
#' fit <- single_release(tags = 803, catch = carphauls$catch,
#'                      recaps = carphauls$tagged,
#'                      prior_recaps = c(108, 181, 82),
#'                      method = "Chapman", unit = "numbers",
#'                      type = 2, tag_mort = 0.3, nat_mort = 0.04,
#'                      chronic_mort = 0.17)
#' summary(fit)
#' 
#' ## Not run: 
#' #######################################
#' # estimate uncertainty with a bootstrap
#' #######################################
#' 
#' boot_fit <- bootstrap(fit, 1000)
#' summary(boot_fit)
#' 
#' ## End(Not run)
single_release <- function(tags, catch, recaps, mean_wt=0, prior_recaps=0, 
                           method, unit, type, tag_mort=0, reporting=1, 
                           nat_mort=0, chronic_shed=0, chronic_mort=0){
  ## check the releases are > recaptures
  if((sum(recaps)+sum(prior_recaps)) > tags) 
    stop("more tagged individuals have been recaptured than were released")
  ## check the method
  if(!method %in% c("Petersen", "Chapman")) 
    stop("incorrect method, currently 'Petersen' and 'Chapman' are implemented")
  ## check units
  if(!unit %in% c("numbers", "kg", "tonnes")) 
    stop("incorrect units must be numbers, kg or tonnes")
  ## so we expect to have more than one
  if(length(catch) != length(recaps)) 
    stop("catch and recaptures must be of the same length")
  ## more checks as they come to mind
  ## if catch and recaps are length 1 we can't bootstrap, 
  ## however this isn't necessarily a problem
  ## define the number of years and check vector lengths
  n_years <- length(prior_recaps)
  ## if some parameters are scalars create vectors with rep
  if(n_years > 1){
    ## generalise with a new function
    if(length(reporting)==1) reporting <- rep(reporting, n_years)
    if(length(nat_mort)==1) nat_mort <- rep(nat_mort, n_years)
    if(length(chronic_shed)==1) chronic_shed <- rep(chronic_shed, n_years)
    if(length(chronic_mort)==1) chronic_mort <- rep(chronic_mort, n_years)
  }
  ## check inputs are correct length, again replace with a function
  if(any(length(reporting)!=n_years, length(nat_mort)!=n_years, 
         length(chronic_shed)!=n_years, length(chronic_mort)!=n_years)) 
    stop("Inputs are not of the same length")  
  ## calculate the adjusted numbers
  ## if there are zero recaptures add an NA
  if(sum(recaps)==0){
    ## if there are no recaptures we don't estimate population size
    N_hat=NA
    var_N=NA
    est<- c(N_hat,var_N)
    names(est)<-c("N_hat","var_N")
   ## apply tag-induced mortality in the first year only
   # assume all releases are lost from the previous season 
    adj_tags <- 0
  }else{
    ## loop over the years to calculate 
    for(i in 1:n_years){
      if(type==1){
        ## calculate the adjusted releases based on tag-induced mortality
        if(i==1){
          ## apply tag-induced mortality in the first year only
          adj_tags <- tags*exp(-tag_mort)
          }
          ## remove recaptures scaled for reporting rate (could improve)
          adj_tags <- adj_tags - (prior_recaps[i] / reporting[i])
          ## apply natural mortality etc
          adj_tags <- adj_tags*exp(-nat_mort[i]-chronic_shed[i]-chronic_mort[i])
      }else if(type==2){
        ## calculate the adjusted releases based on tag-induced mortality
        if(i==1){
          ## apply tag-induced mortality in the first year only
          adj_tags <- tags*exp(-tag_mort)
        }
        ## remove half recaptures scaled for reporting rate (could improve)
        adj_tags <- adj_tags - 0.5*(prior_recaps[i] / reporting[i])
        adj_tags <- adj_tags*exp(0.5*(-nat_mort[i]-chronic_shed[i]-chronic_mort[i]))
        adj_tags <- adj_tags - 0.5*(prior_recaps[i] / reporting[i])
        adj_tags <- adj_tags*exp(0.5*(-nat_mort[i]-chronic_shed[i]-chronic_mort[i]))
      }else stop("type of fishery is not defined")
    ## adjust recaptures in the survey season by dividing by reporting rate 
    adj_recaps <- sum(recaps) / reporting[n_years]
    ## now check that the adjusted releases are > 0
    if(adj_tags <= 0){
      ## alternately N_hat could be zero
      N_hat <- NA 
      var_N <- NA
      est<-c(N_hat,var_N)
      names(est)<-c("N_hat","var_N")
      warning("predicted releases are less than zero based on other parameters")
    }else if(method=="Petersen"){
      ## Petersen numbers and weight are the same
      est <- petersen(adj_tags, sum(catch), adj_recaps)
    }else if(method=="Chapman" & unit %in% c("numbers")){
      ## Chapman numbers
      est <- chapman_n(adj_tags, sum(catch), adj_recaps)
    }else if(method=="Chapman" & unit %in% c("kg", "tonnes")){
      ## Chapman weight
      est <- chapman_wt(adj_tags, sum(catch), adj_recaps, mean_wt)
    }else stop("method and unit combination not available")
    }
  }
  ## store the data
  ###*******************************************###
  hauls <- data.frame(cbind(catch, recaps))
  ## wrap the results in a list
  obj <- list("Estimate" = est,
              "Method" = method,
              "Unit" = unit,
              "Type" = type,
              "Releases" = tags,
              "Hauls" = hauls,
              "MeanWeight" = mean_wt,
              "PriorRecaps" = prior_recaps,
              "TagMort" = tag_mort, 
              "Reporting" = reporting,
              "NatMort" = nat_mort,
              "ChronicShed" = chronic_shed,
              "ChronicMort" = chronic_mort,
              "TagsAvailable" = adj_tags)
  ## add class 'srelease'
  class(obj) <- 'srelease'
  obj
}

## ************************************************ ##
## S3 Generics

#' @export
#' @rdname bootstrap
bootstrap.srelease <- function(x, nboot, ...){
  ## check the there are sufficient rows in the data
  if(nrow(x$Hauls) <= 1) 
    stop("must be more than one haul to undertake bootstrap")
  if(nrow(x$Hauls) < 10) 
    warning("there are less than 10 hauls used in the bootstrap")
  ## extract the various components of x
  tags <- x$Releases
  hauls <- x$Hauls
  prior_recaps <- x$PriorRecaps 
  tag_mort <- x$TagMort
  reporting <- x$Reporting
  nat_mort <- x$NatMort
  chronic_shed <- x$ChronicShed
  chronic_mort <- x$ChronicMort
  ## extract the number of years to loop over
  n_years <- length(prior_recaps)
  ## object to store the bootstrapped estimates
  boot_est <- data.frame(matrix(NA, nrow=nboot, ncol=3))
  names(boot_est) <- c("boot_est", "boot_var", "boot_recaps")
  ## loop over the number of simulations 
  for(j in 1:nboot){
    ## stochastically apply tag loss and mortality by season
    for(k in 1:n_years){
      ## in the first year we apply tag induced mortality
      if(k==1){
        adj_tags <- rbinom(n=1, size=tags, prob=(exp(-tag_mort)))
      }
      if(x$Type == 1){
        ## remove the recaptures scaled for reporting rate and rounding up
        adj_tags <- ceiling(adj_tags - (prior_recaps[k] / reporting[k]))
        ## apply tag shedding and natural mortality 
        adj_tags <- rbinom(n=1, size=adj_tags, 
                           prob=exp(-nat_mort[k]-chronic_shed[k]-chronic_mort[k]))
      }else if(x$Type == 2){
        ## remove half recaptures scaled for reporting rate and rounding up
        adj_tags <- ceiling(adj_tags - 0.5*(prior_recaps[k] / reporting[k]))
        ## apply half tag shedding and natural mortality 
        adj_tags <- rbinom(n=1, size=adj_tags, 
                           prob=exp(0.5*(-nat_mort[k]-chronic_shed[k]-chronic_mort[k])))
        ## remove second half of the recaptures
        adj_tags <- ceiling(adj_tags - 0.5*(prior_recaps[k] / reporting[k]))
        ## apply half tag shedding and natural mortality 
        adj_tags <- rbinom(n=1, size=adj_tags, 
                           prob=exp(0.5*(-nat_mort[k]-chronic_shed[k]-chronic_mort[k])))
      }else stop("type of fishery is not defined")
    }
    ## check adj_tags >=1 otherwise don't estimate population size and return NA
    if(adj_tags <1){
      boot_est[j,] <- c(NA, NA, NA)
    }else{
      ## resample the haul data 
      j_sample <- sample(nrow(hauls), replace=TRUE)
      j_recaps <- sum(hauls[j_sample,]$recaps)/ reporting[n_years]
      j_catch  <- sum(hauls[j_sample,]$catch)
      ## calculate Chapman estimate
      if(x$Method=="Petersen"){
        ## Petersen numbers
        boot_est[j,] <- c(suppressWarnings(petersen(adj_tags, j_catch, j_recaps)), j_recaps)
      }else if(x$Method=="Chapman" & x$Unit %in% c("numbers")){
        ## Chapman numbers
        boot_est[j,] <- c(suppressWarnings(chapman_n(adj_tags, j_catch, j_recaps)), j_recaps)
      }else if(x$Method=="Chapman" & x$Unit %in% c("kg", "tonnes")){
        ## Chapman weight
        boot_est[j,] <- c(suppressWarnings(chapman_wt(adj_tags, j_catch, j_recaps, x$MeanWeight)), j_recaps)
      }
    }
  }
  ## store the inputs and results 
  obj <- list("srelease_obj" = x,
              "Boot_estimates" = boot_est)
  ## add a class
  class(obj) <- "bsamples"
  ## return the results
  return(obj)
}

#' Print method for single release, single recapture models
#'
#' S3 method for single release, single recapture models 
#' tag recapture models of abundance
#' @param x object of class srelease
#' @param ... additional parameters
#' @export
print.srelease <- function(x, ...){
  ## define the model used
  cat(x$Method, "estimate of population", x$Unit, "\n")
  print(x$Estimate)
}

#' Summary method for single release, single recapture models
#'
#' S3 method for single release, single recapture models 
#' tag recapture models of abundance
#' @param object object of class srelease
#' @param ... additional parameters
#' @export
summary.srelease <- function(object, ...){
  x <- object
  ## print the model inputs 
  cat(x$Method, "estimate of", x$Unit, "\n")
  print(x$Estimate)
  cat("With cv", sqrt(x$Estimate["var_N"]) / x$Estimate["N_hat"], "\n")
  ## summarise the inputs
  cat(x$Releases, "tags were released and", sum(x$PriorRecaps),
      "were subsequently recaptured prior to this season \n")
  cat(sum(x$Hauls["catch"]), x$Unit, "were captured in the current survey and",
      sum(x$Hauls["recaps"]), "were tagged \n")
  cat("The following parameters were specified \n")
  cat("Initial tag-induced mortality =", x$TagMort, "\n")
  cat("Tag reporting rates by season =", x$Reporting, "\n")
  cat("Natural mortality by season =", x$NatMort, "\n")
  cat("Chronic tag shedding by season =", x$ChronicShed, "\n")
  cat("Chronic tag-induced mortality by season =", x$ChronicMort, "\n")
  cat("Tags available by season =", x$TagsAvailable,"\n")
}

#' S3 method for bootstrapped confidence intervals
#'
#' Calculate bootstrapped confidence intervals for object of
#' class bsamples
#' @param object object of class bsamples
#' @param quantiles bootstrap quantiles (default 0.025, 0.5, 0.975)
#' @param ... additional parameters
#' @export
summary.bsamples <- function(object, quantiles=c(0.025, 0.5, 0.975), ...){
  ## define the quantiles
  quants <- quantile(object$Boot_estimates$boot_est, probs=quantiles, na.rm=TRUE)
  se <- sd(object$Boot_estimates$boot_est, na.rm=TRUE)
  cv <- se / object$srelease_obj$Estimate["N_hat"]
  names(se) <- "boot_SE"
  names(cv) <- "boot_CV"
  ## extract the information regarding the 
  ## print the model inputs 
  cat(object$srelease_obj$Method, "estimate of", object$srelease_obj$Unit, "\n")
  cat(object$srelease_obj$Releases, "tags were released and", 
      sum(object$srlease_obj$PriorRecaps), 
      "were subsequently recaptured in previous seasons \n")
  cat(sum(object$srelease_obj$Hauls$catch), object$srelease_obj$Unit,
      "were captured in the current survey and ", 
      sum(object$srelease_obj$Hauls$recaps), " were tagged \n")
  cat("The following parameters were specified \n")
  cat("Initial tag-induced mortality =", object$srelease_obj$TagMort, "\n")
  cat("Tag reporting rates by season =", object$srelease_obj$Reporting, "\n")
  cat("Natural mortality by season =", object$srelease_obj$NatMort, "\n")
  cat("Chronic tag shedding by season =", object$srelease_obj$ChronicShed, "\n")
  cat("Chronic tag-induced mortality by season =", 
      object$srelease_obj$ChronicMort, "\n")
  ## construct the output
  out <- c(object$srelease_obj$Est[1], se, cv, quants)
  #names(out) <- c("Estimate", "lower", "upper")
  out
}
