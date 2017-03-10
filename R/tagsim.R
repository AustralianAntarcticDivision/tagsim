#' tagsim: A simple population model for tag-based simulations
#'
#' A simple population model for tag-based simulations for fisheries based on
#' the fisheries operating model written by Philippe Ziegler.
#'
#' Simulation model
#'
#' The simulation model is aimed at evaluating tag-based scenarios, however
#' scenarios without tagging data can be implemented.
#'
#' Assessment model
#'
#' Simple assessments are implemented including
#'
#' 1) Tag-based estimates
#'
#' Simulation Processes
#'
#' Processes in the model occur in the following order
#'
#' 1 New season processes
#'
#' 1.1 create a temp vector to store the population and tag numbers
#' last year
#'
#' 1.2.1 Move untagged population
#'
#' 1.2.2 Move the tagged population
#'
#' 1.3 Estimate recruitment from last season numbers (there is no growth)
#'
#' 1.3.1 add the recrutiment (recruits don't move) consider saving this
#'
#' 2 Harvest and natural mortality
#'
#' 2.1 calculate population size by strata (tagged + untagged)
#'
#' 2.1.1 update the harvest rate/strategy based on last seasons abundance
#'
#' 2.1.2 calculate an F based on fishery type (Ricker I or II)
#'
#' 2.2 calculate the total catch & natural mortality
#'
#' 2.3 remove from the population tagged and untagged
#'
#' 2.4 calculate the number of tags observed etc
#'
#' 3 Assign the end of season population size to the storage
#'
#' 4 Estimate abundance and store the result
#'
#' 5 move to the next year and repeat
#'
#' @docType package
#' @name tagsim
NULL
#> NULL
