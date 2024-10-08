% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_study_opt_ar1.R
\name{sim_study_opt_ar1}
\alias{sim_study_opt_ar1}
\title{Simulation study using optimization routine, multiple cores}
\usage{
sim_study_opt_ar1(
  outDir = here::here("out.txt"),
  betaVec = c(0, 0.12, 0.25, 0.5),
  rhoVec = c(0, 0.25, 0.5, 0.75, 0.95),
  sigmaVec = c(0.25, 0.5, 0.75),
  nTVec = c(10, 25, 50),
  nSims = 200,
  nBootSims = 500,
  setSeed = NULL,
  nCores = NULL,
  missing = F
)
}
\arguments{
\item{outDir}{Character string. Path to output file}

\item{betaVec}{Numeric vector. Values for beta_1 (slope/trend parameter)}

\item{rhoVec}{Numeric vector. Values for autoregressive parameter}

\item{sigmaVec}{Numeric vector. Values of standard deviation of noise}

\item{nSims}{Numeric scalar. Number of time series to simulate}

\item{nBootSims}{Numeric scalar. Number of bootstrap data sets}

\item{setSeed}{Numeric scalar. Value of the seed for simulations. (Default = NULL, a random number between 1-e7 is selected)}

\item{nCores}{Numeric scalar. Specify the number of cores to utilize (Default = NULL, utilizes n-1 cores)}

\item{missing}{Boolean. Whether to simulate with missing data (Default = F).
If T then a single missing value is added at random to the response}

\item{nVec}{Numeric vector. Values for the length of time series to simulate}
}
\description{
Performs a simulation study to assess the performance of the test.
Utilizes multiple cores to spread bootstrap samples over multiple cores
}
\details{
For convenience the intercept, beta_0 is set to zero
}
\examples{
\dontrun{
}


}
