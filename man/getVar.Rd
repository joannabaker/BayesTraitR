% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manipulateLogs.R
\name{getVar}
\alias{getVar}
\title{Extract variance from BayesTraits output.}
\usage{
getVar(lf)
}
\arguments{
\item{lf}{The path and filename of the log file the variance is wanted for.}
}
\value{
List with elements:
\item{BG}{a vector of the extracted phylogenetic variance}
\item{dBG}{computed density distribution}
\item{dLogBG}{computed logged density distribution}
}
\description{
This will work with raw output.
}
