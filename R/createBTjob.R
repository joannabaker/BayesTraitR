#' @title Generate DistData file
#' @description	Takes a data file and generates a DistData file for
#'     analysis in BayesTraits with multiple datapoints for some tips.
#' @param unlinked The data table that is wanted for further analysis in BayesTraits.
#'     This MUST contain taxa names or row IDs in the first column.
#' @param link Currently, the value must be "Linked" or "Unlinked".
#'      This will be applied across all taxa. Future updates may allow
#'      variable values across taxa.
#' @return An R data frame that can be exported to a text file and used as input
#' for the DistData command in BayesTraits.This function is used internally by
#' [createBTjob()].
#' @export

linkvalues = function(unlinked, link = "Linked"){
  name = colnames(unlinked)[1]
  taxa = unique(unlinked[,name])

  cols = colnames(unlinked)
  cols = cols [cols != name]

  linked = data.frame(matrix(nrow = 0, ncol = ncol(unlinked)))
  colnames(linked)= colnames(unlinked)
  #linked = unlinked[match(unique(unlinked[,name]), unlinked[,name]),]
  for(tx in taxa){
    tmp = unlinked[unlinked[,name] == tx,]
    if(nrow(tmp) == 1) next
    new = tmp[1,]
    for( i in 1:length(cols)){
      vals = tmp[,cols[i]]
      new[,cols[i]] = paste(vals[!is.na(vals)], collapse = ",")
    }

    linked = rbind(linked, new)
  }

  linked$linked = link
  linked= linked[,c(1,ncol(linked), 2:(ncol(linked)-1))]
  return(linked)

}

#' @title Generate priors for Beta estimation in BayesTraits.
#' @description Generate priors formatted for use as an 'optional argument' in [createBTjob()].
#'     This is specifically designed to generate priors for regression betas in order to
#'     avoid typing out variable names manually.
#' @param data The data table that is wanted for further analysis in BayesTraits.
#'     This must contain ONLY the columns required for subsequent analysis,
#'     and must have names matching the tree as the first column.
#' @param prior Name of the prior to apply. Can only be one of "normal", "unif",
#'     "gamma", "exp". See BayesTraits manual for details.
#' @param pars Numeric values defining the limits/parameters of the prior.
#'     See BayesTraits manual for details.
#' @return A character string used in BayesTraits input files to define priors
#'     for a set of Betas. Designed for use with [createBTjob()].
#' @export
BTBetapriors = function(data, prior, pars) {
  priors = NULL;
  for(i in 3:ncol(data))
    priors = c(priors, paste0("prior Beta-", i-2, " ", prior, " ", paste(pars, collapse = " ")))
  return (priors)}


#' @title Generate BayesTraits analysis files.
#' @description Flexible function for automatically generating all of the files
#'     required for running phylogenetic statistical comparative analyses in BayesTraits.
#' @param fm An object of class \code{\link[stats]{formula}} that is a symbolic description of the
#'     desired model to be fitted in BayesTraits. See details.
#' @param dataset The data table that is wanted for further analysis in BayesTraits.
#'     This must contain ONLY the columns required for subsequent analysis,
#'    and must have names matching the tree as the first column.
#' @param tree A single tree of class \code{"phylo"} containing tip-names as found in the first
#'    column of \code{"dataset"}.
#' @param jobname A unique character string used for identifying the model.
#'     All created files will have this tag.
#' @param bi Numeric value determining the amount of iterations to remove from
#'     the beginning of the chain as burn-in. Defaults to \code{100000}.
#' @param it Numeric value determining the total number of iterations to run.
#'     Default is \code{1100000}.
#' @param sa Numeric value determining the "thinning" or "sampling" parameter -
#'     i.e. how many iterations to run between each sample. Default is \code{1000}.
#' @param model A required argument that specifies which model is to be run in
#'     BayesTraits. See \href{http://www.evolution.reading.ac.uk/BayesTraitsV4.0.1/BayesTraitsV4.0.1.html}{BayesTraits} manual for list of available models.
#' @param MCMC Boolean operator that determines whether the analysis to be run
#'     will be run using Markov chain Monte Carlo methods (\code{MCMC = T}) or in Maximum Likelihod (\code{MCMC=F}).
#'     Defaults to \code{MCMC=TRUE}.
#' @param reps Numeric value (1:999) defining the number of replicates of this model to be run.
#'     Defaults to 1.
#' @param optarg an optional character vector with arguments to
#'     include in the input file for running BayesTraits. This includes commands
#'     like \code{"varrates"} or \code{"lambda"} in addition to prior specifications and
#'     stepping stone sampling. See \href{http://www.evolution.reading.ac.uk/BayesTraitsV4.0.1/BayesTraitsV4.0.1.html}{BayesTraits} manual for list of available commands.
#' @param outdir Path to desired folder where all generated job files will be saved.
#'     If this parameter is left unspecified, files will be generated in the current
#'     working directory.
#' @param DistData Boolean operator that specifies whether there are sampled
#'     values to be accounted for (\code{DistData=TRUE}), or if each tip has only a single
#'     corresponding data point (\code{DistData=FALSE}). This allows the user to input a sample of data
#'     for a taxon/range of taxa instead of a single point. This requires the input
#'     data file to contain multiple entries for at least one of the sampled taxa.
#' @param link If \code{DistData=TRUE}, the user must specify whether samples are "linked"
#'     (traits are taken from the same source or individual) or "unlinked" (samples)
#'     are not taken from a single source or individual.
#' @param contrasts A list, whose entries are values (numeric matrices, \link[stats]{formula}s,
#'  or character strings naming functions) to be used as replacement values for the \link[stats]{contrasts}
#'  replacement function and whose names are the names of columns of \code{dataset} containing \link[base]{factor}s.
#'  Argument and description from \link[stats]{model.matrix}.
#'  @param names.col An argument specifying the column of the dataset that links to the names in the tree.
#'  If left unspecified, will search for names in the first column.
#' @return Generates files in the current working directory that can be used
#'     to run BayesTraits analyses.
#' @details BayesTraits is in constant development and as such not all eventualities are
#'     accounted for in this package. This function is designed with the informed
#'     user in mind. It is assumed that the user knows which model and parameter
#'     combinations are appropriate for their own data.
#'
#'     Using this function, models can be specified symbolically. The response
#'     variable must always be a numeric vector. Models can be specified as described
#'     in the documentation for [lm()].
#'
#' @export
#' @importFrom ape drop.tip write.nexus
#' @importFrom utils read.table write.table
#' @importFrom stats model.matrix
createBTjob <- function(fm, dataset, tree, jobname = "BTjob", bi = 100000, it = 1100000,
                        sa = 1000, model, MCMC = T, reps = 1, optarg, outdir, DistData, link,
                        contrast_arg = NULL, names.col = colnames(dataset)[1]){
  # Specify output directory
  if(missing(outdir)) outdir = "."

  # Quick check that the names of the species are in the tree
  if(!any(dataset[,names.col] %in% tree$tip.label))
    stop("No species in dataset found in tree. Ensure names.col is specified and that it includes names corresponding to tip labels.")

  # For contrast coding to work as expected, categorical variables should be treated as factors.
  vartype = sapply(all.vars(fm),function(x)class(dataset[,x]))
  if(any(vartype == "character")){
    warning("Character variables detected; converting to factors.")
    dataset[,names(vartype)[vartype == "character"]] <- lapply(dataset[names(vartype)[vartype == "character"]] , factor)
    vartype = sapply(all.vars(fm),function(x)class(dataset[,x]))
  }

  # Create and write the tree file
  .mistx = tree$tip.label[!tree$tip.label %in% dataset[,names.col]]
  if(length(.mistx) > 0){
    cat("Removing", length(.mistx), "tips from the phylogeny:\n",
        paste0(.mistx, collapse = "\n "), "\n\n")
    tree = drop.tip(tree, .mistx)
  }
  write.nexus(tree, file = paste0(outdir, "/",jobname, ".trees"))

  # Ensure all taxa are in the data file
  .mistx = dataset[,1][!dataset[,names.col] %in% tree$tip.label]
  if(length(.mistx) > 0){
    cat("Removing", length(.mistx), "rows from the dataset:\n",
        paste0(.mistx, collapse = "\n "), "\n")
    dataset = dataset[!dataset[,names.col] %in% .mistx,]
  }

  # Calculate the model matrix
  mf = model.matrix(fm, dataset, na.action = 'na.pass', contrasts.arg = contrast_arg)

  # Identify the Y variable and modify the design matrix to BT input format
  mf[,1] = dataset[,all.vars(fm)[1]]

  # Add the names column to the input data
  modeldata = cbind.data.frame(dataset[,names.col],mf)

  # Modify the column names of the dataset for the intercept / names column
  colnames(modeldata)[1:2] = c(names.col,all.vars(fm)[1])

  ## I have commented this out as it's not really good practice
  ## But I found it helpful for interpreting the outputs of sum contrasts vs dummy codes
  # # Modify the column names for any factors
  # if(any(vartype == "factor")){
  #   factors = vartype[vartype == "factor"]
  #   for(f in 1:length(factors)){
  #     fac = names(factors)[f]
  #     faclvls = levels(dataset[,fac])
  #
  #     # Remove group name from column
  #     colnames(modeldata) = gsub(fac, "", colnames(modeldata))
  #
  #     # If we had per-group coding, names should be OK
  #     if(any(faclvls %in% colnames(modeldata))) next
  #
  #     # If we have contrast coding like sum diffs, adjust name
  #     for(L in 1:length(faclvls)){
  #       if(L %in% colnames(modeldata))
  #         colnames(modeldata)[which(colnames(modeldata) == L)] = as.character(faclvls)[L]}
  #
  #
  #   }}

  # BayesTraits does not accept spaces or special characters in column names
  colnames(modeldata) = gsub("-|//*| |:", "_", colnames(modeldata))

  # Chain settings
  if(MCMC==F) conditions ="" else
    conditions = c(paste("burnin", format(bi, scientific =F)),
                   paste("iterations", format(it, scientific = F)),
                   paste("sample", format(sa, scientific = F)))

  # Add optional arguments
  if(missing(optarg)) optarg=""

  # Add distribution data table if required
  if(!missing(DistData)){

    # If we have no duplicates, this is a mistake. Check with user.
    if(nrow(dataset) == length(unique(dataset[,1])))
      stop("DistData specified with no sampled values found in dataset.")

    # Secondly, check link specification
    if(missing(link)){
      warning("DistData specified with no link value defined.
                               Defaulting to Linked.")
      link = "Linked"}

    # For any samples of data, create a DistData table
    linkeddata = linkvalues(dataset,link)
    dd = paste0(jobname, "-DistData.txt")
    write.table(linkeddata, file = paste0(outdir, "/", dd), sep  = "\t",
                col.names = F, row.names = F, quote = F)

    # Create a new input file (we need only one data point per taxon)
    dataset = dataset[match(unique(dataset[,1]), dataset[,1]),]

    # Add command to optional arguments
    optarg = c(optarg, paste0("DistData ", dd))

  }

  # Create replicates and save them to the output folder
  for(i in 1:reps){
<<<<<<< HEAD
    write.table(modeldata, file = paste0(outdir, "/",jobname, "-", stringr::str_pad(i, 3, pad =0), ".txt"),
=======
<<<<<<< HEAD

=======
<<<<<<< HEAD
=======
    print(i)
>>>>>>> de3b668c7b81c478b5983f284051de7600d7b5eb
>>>>>>> 9a49a23ea4a4a01e39978b5a70253494892fba05
    write.table(modeldata, file = paste0(jobname, "-", stringr::str_pad(i, 3, pad =0), ".txt"),
>>>>>>> bc64730fcc96136cd2ca19819ca657ba6f2d1e2a
                sep = "\t", col.names = T, row.names = F, quote = F)}

  # Create and save the input file
  inf = c(model, ifelse(MCMC==T,2,1), optarg,conditions,"run")
  writeLines(inf, con = paste0(outdir, "/", jobname, ".infile"))
}

