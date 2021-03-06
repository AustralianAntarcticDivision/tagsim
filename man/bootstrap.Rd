% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_tagr.R, R/multi_release.R,
%   R/single_release.R
\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap.srelease}
\alias{bootstrap.mrelease}
\alias{bootstrap.srelease}
\title{Bootstrap method for single release tag-return studies}
\usage{
bootstrap(x, nboot, ...)

\method{bootstrap}{mrelease}(x, nboot, ...)

\method{bootstrap}{srelease}(x, nboot, ...)
}
\arguments{
\item{x}{an object of class srelease}

\item{nboot}{number of bootstrap samples}

\item{...}{additional parameters}
}
\description{
Bootstrap method for single release tag-return studies
}
\details{
Estimate confidence intervals using a non-parametric bootstrap. This method
incoporates uncertainty in tag-induced mortality, natural mortality and tag 
shedding as a series of Bernoulli trials. @seealso \code{\link{single_release}}
for more details.
}
