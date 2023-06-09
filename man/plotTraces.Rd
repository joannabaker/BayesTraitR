% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manipulateLogs.R
\name{plotTraces}
\alias{plotTraces}
\title{Plot parameter values against index value.}
\usage{
plotTraces(files, cols = "all", table = T)
}
\arguments{
\item{files}{path or filename of the file a summary is required for.
This can be a list of filenames - if so, all identical values in each
table will be plotted on a single chart. This is useful for comparing estimates
of the same parameter across multiple replicate,s for example.
File must be a tab-delimited text file with column headers.
This function is designed to be compatible with the output of \code{\link[=trimmedlogs]{trimmedlogs()}}.}

\item{cols}{A string or list of strings defining the column names or indices
for which summary information is required. If unspecified, takes the
value "all" which will return traces for the following columns
(if present): Likelihood, Intercept, Slope(s)/Beta parameters,
Variance, R-suqred, Number of local scalars (branch or node), and Lambda.}

\item{table}{if TRUE, then the function will accept a list of \code{data.frames}
rather than a list of filenames.}
}
\description{
Useful for visualizing parameter traces from MCMC analyses in
BayesTraits but can in principle be applied to any data.
}
