% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppBTLogs.R
\name{predBTlog}
\alias{predBTlog}
\title{Calculate predictions from BayesTraits regressions.}
\usage{
predBTlog(input, output = NULL, plot = F, burnin = 0)
}
\arguments{
\item{input}{A character vector defining the original input file (must have column headers).
Note that this can be a modified version of the original input e.g. if one wants to fix a value at a mean - as long as columns are retained in identical order and format.}

\item{output}{an optional argument. If defined, a character vector locating the original output file (unmodified).
If unspecified, the function will search for the original input file as if it had been run through BayesTraits, appending .Log.txt.}

\item{burnin}{A value specifying the number of rows to remove from the beginning of each log file. Default value is 0 and so the entire chain will be used.}
}
\value{
A matrix of predicted values for terminal taxa (if original input was MCMC, one iteration per column)
}
\description{
Function to calculate predictions from BayesTraits input/output files. Incompatible with other programs.
}
