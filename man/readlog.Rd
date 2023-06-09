% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manipulateLogs.R
\name{readlog}
\alias{readlog}
\title{Read in BayesTraits Log Files (raw)}
\usage{
readlog(file)
}
\arguments{
\item{file}{Takes a single string defining the direct path to a single log file.}
}
\value{
R data.frame object that contain the trimmed log file.
}
\description{
Removes the header from BayesTraits output files.
Can only work with log files (.txt.Log.txt) as directly output from BayesTraits models 1-4, 7 and 9.
Will not work with output from any other program or model.
}
