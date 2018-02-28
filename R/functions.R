#' Baranov catch equation
#'
#' Returns catch for given F for population size N
#' with natural and fishing mortality M and F respectively
#'
#' Based on code written by Bill de la Mare
#' @param N population size
#' @param nat_mort instantaneous natural mortality
#' @param fish_mort instantaneous fishing mortality
#' @export
baranov_catch <- function(N, nat_mort, fish_mort)
{
  N*fish_mort/(nat_mort+fish_mort)*(1-exp(-nat_mort-fish_mort))
}

#' Calculates catch for a given F
#'
#' Returns catch from N for a given F, used in \code{catchtoFfinding} the F that
#' corresponds to a given catch.
#'
#' Based on code written by Bill de la Mare
#' @param N population size
#' @param nat_mort instantaneous natural mortality
#' @param fish_mort instantaneous fishing mortality
#' @export
F_to_catch <- function(N, nat_mort, fish_mort){
  ## calculate total mortality Z
  Z <- fish_mort + nat_mort
  # calculate the fraction of the population that dies from fishing
  N * (1-exp(-Z))*fish_mort/Z
}

#' Calculate F for a given catch
#'
#' Root of this function is the F that corresponds to a given catch
#'
#' Based on code written by Bill de la Mare
#' @param N population size
#' @param nat_mort instantaneous natural mortality
#' @param fish_mort instantaneous fishing mortality
#' @param catch catch
#' @export
catch_to_F <- function(N, nat_mort, fish_mort, catch)
{
  catch - F_to_catch(N, nat_mort, fish_mort)
}

