% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_release.R
\name{pars}
\alias{pars}
\title{Parameters for specifying multiple release tag model}
\arguments{
\item{mean_wt}{mean weight of a single fish (optional argument for Chapman 
weight method)}

\item{method}{The method to use for calculating population size (options are
"Petersen" or "Chapman")}

\item{unit}{the unit of measurement (either "numbers", "kg" or "tonnes")}

\item{type}{does fishing occur over a short period (Ricker type 1) or is
extended over the year and competes with natural mortality and tag sheddding
(Ricker type 2)}

\item{tag_mort}{vector of initial tag-induced mortality}

\item{reporting}{vector of annual tag reporting rates}

\item{nat_mort}{vector of natural mortality (instantaneous)}

\item{chronic_shed}{vector of chronic (ongoing) tag shedding}

\item{chronic_mort}{vector of chronic (ongoing) tag-induced mortality}
}
\description{
Parameters for specifying multiple release tag model
}
\details{
A list that specifies the following parameters. Mixture of text, integer and
vectors
}
\seealso{
\code{\link{multi_release}} for more details.
}
