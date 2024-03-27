#' @title Generate BayesPhylogenies analysis files.
#' @description Flexible function for automatically generating the BayesPhylogenies
#'     command block required for running a single-partitioned model.
#' @param gamma Add gamma rate heterogeneity with discrete classes = \code{gamma}.
#' @param chainsN Specify the number of chains to be run. Default value = \code{1}.
#' @param model A character string specifying the model of evolution. It must be one of
#'     "GTR" (the default), "GTNR", "HKY85", "F81", "SYM", "JC", "K2P" for sequence data.
#'     For binary data, model can be "M1P", "M2P" and for multi-state data the "KSTATES"
#'     model can be used.
#' @param burnin Numeric value determining the amount of iterations to remove from
#'     the beginning of the chain as burn-in. Defaults to \code{100}.
#' @param iterations Numeric value determining the total number of iterations to run.
#'     Default is \code{1100}.
#' @param sample Numeric value determining the "thinning" or "sampling" parameter -
#'     i.e. how many iterations to run between each sample. Default is \code{500}.
#' @param seed A numeric value specifying the random number seed for a model.
#'     If left unspecified, BayesPhylogenies will randomly specify this. If specified,
#'     it is possible to replicate a Markov chain using a given seed.
#' @param autorun A Boolean operator. If \code{autorun=TRUE}, then the program
#'     will run immediately without user-input if the output of this function
#'     is executed within BayesPhylogenies.
#' @param rjpatterns Boolean operator specifying whether or not to implement the
#'     Pagel and Meade (2004) reversible-jump mixture model for pattern heterogeneity.
#'     Defaults to \code{FALSE}.
#' @param blsets An optional numeric value that will implement the heterotachy model
#'     as described in (Meade and Pagel, 2008; Pagel and Meade, 2008) for a given
#'     number of branch length sets. Incompatible with \code{rjbranchlengths}.
#' @param rjbranchlengths an optional numeric value that implements the reversible-jump
#'     heterotachy model for a given number of branch length sets. Incompatible with
#'     \code{blsets}.
#' @param topologies An optional numeric value that allows for multiple topologies
#'     to be estimated. If unspecified, the standard single topology model is run.
#'     It is not recommended that this value be larger than 2.
#' @param optionalarguments An optional character vector containing further arguments
#'     to be passed to BayesPhylogenies. Please see program documentation for options
#'     and examples of accepted arguemnts.
#' @return A character string representing the input to BayesTraits
#' @details BayesPhylogenies is in constant development and as such not all eventualities are
#'     accounted for in this package. This function is designed with the informed
#'     user in mind. It is assumed that the user knows which model and parameter
#'     combinations are appropriate for their own data.
#'
#'
#' @export
BP_SPblock = function(chainsN=1, model = "GTR", iterations = 1100, sample = 500,
                         burnin = 100, seed = "random",
                         autorun=T, rjpatterns=F, gamma,
                         blsets,rjbranchlengths, topologies, optionalarguments){

  # Construct template BayesPhylogenies block for a single-partitioned model
  init = paste0("Begin BayesPhylogenies;\n\t#Create partition\n\t\tChains ", chainsN,";\n\t\tCreatePart\t Default", "\t", model,";")

  # Add gamma
  if(!missing(gamma)) init = c(init, paste("\t\tgamma", gamma))

  # Construct MCMC chain settings
  init = c(init, paste("\n\n\t#Set up MCMC chain\n\titerations",
                       format(iterations, scientific =F), "\n\tsample",
                       format(sample, scientific = F), "\n\tburnin",
                       format(burnin, scientific=F)))

  if(seed !="random") init = c(init, paste("\tseed", seed))

  # RJ pattern heterogeneity?
  if(rjpatterns) init = c(init, "\n\t#Set up RJ pattern heterogeneity\n\trjpatterns")

  # Heterotachy
  if(!missing(rjbranchlengths) & !missing(blsets)) stop("Multiple heterotachy model types specified.") else{
    if(!missing(blsets))init = c(init, paste("\n\t#Set up heterotachy (not RJ)\n\tblsets", blsets))
    if(!missing(rjbranchlengths))init = c(init, paste("\n\t#Set up RJ heterotachy \n\trjbranchlengths", rjbranchlengths))
  }

  # Multiple topologies
  if(!missing(topologies)) init = c(init, paste("\n\t#Set up RJ topologies\n\ttopologies", topologies))

  # Optional arguments
  if(!missing(optionalarguments)) init = c(init, paste0("\n\t", optionalarguments, collapse = "\t"))

  #autorun
  if(autorun) init =c(init, "\n\t#Run model\n\tautorun")

  # end
  init = c(init, "\nEnd;")

  # Output
  return(init)
}
