## TODO: Account for equaltrees analyses in post-processing tools

#' @title Read in BayesTraits Log Files (raw)
#' @description Removes the header from BayesTraits output files.
#'     Can only work with log files (.txt.Log.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not work with output from any other program or model.
#' @param file Takes a single string defining the direct path to a single log file.
#' @return R data.frame object that contain the trimmed log file.
#' @export

readBTlog = function(file){
  skip = (which(grepl("Tree No", readLines(file)))) - 1
  log = read.table(file, sep = "\t", header = T, stringsAsFactors = F,
                   comment.char = "*", skip = skip, fill = T)
  return(log)
}

## TODO: Fix file specification and parameter input for BTlogs

#' @title Trim header from BayesTraits log files.
#' @description Removes the header from BayesTraits output files.
#'     Can only work with log files (.txt.Log.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not work with output from any other program or model.
#' @param dir Directory of files for which we want trimmed output. Defaults to the current working directory.
#' @param out The name of the folder to which we want output saved to. Defaults to the current working directory.
#' @param pat If no files are specified, the function will use the string defined by pat to search for files in dir. If only one or a few specific files are desired, then see the \code{files} parameter.
#' @param start A unique search string that defines the start of the columns of output (header row).
#'     Defaults to Tree No which should work for the models mentioned above - check if testing other BayesTraits outputs.
#' @param tail Takes a single numeric value.
#'     If defined, tail defines the number of iterations to be retained from the END of the output log.
#'     Defaults to NA, which will return all rows in the original log file (though see the resample parameter)
#'     If tail = 1000, the output will return only the last 1000 rows.
#' @param resample Takes a single numeric value.
#'     If defined, the function will resample this number of rows from the original log file.
#'     Note that this number must be greater than the total number of rows in the original file, or this will fail.
#' @param files Takes a string or list of strings defining specific log files to be trimmed.
#'     If this is defined, \code{dir} and \code{pat} are ignored.
#' @return A \code{list} containing the trimmed log file, each as a \code{data.frame}.
#' @export


trimBTlog = function(dir = ".", out = ".", pat = ".txt.Log.txt", start = "\tTree No", tail = NA, resample = NA, files= NULL){

  # Set working directory
  home = getwd()
  setwd(dir)

  # Detect all log files
  if(is.null(files)){
    logs = list.files(pattern = pat)
  } else logs = files


  # Generate folder for trimmed log files
  dir.create(out)

  # Generate object for returned log data frames
  results = list ()

  # Loop through and strip the header material out of each log file.
  for(l in 1: length(logs))
  {

    # Read in log file
    .tmp = readBTlog(logs[l])

    # remove NA columns for clarity
    .tmp = .tmp[,!apply(.tmp,2,function(x)all(is.na(x)))]
    .tmp = .tmp[,!apply(.tmp,2,function(x)all(x == "N\\A"))]

    if(!is.na(tail)){
      if(nrow(.tmp) < tail) stop("Number of requested samples is higher than the total number of resampled entries. Check requested values again.")
      .tmp = .tmp[(nrow(.tmp) - tail+1):nrow(.tmp),]
    }

    # SaAve the trimmed output file
    if(!is.na(resample)) {
      if(resample > nrow(.tmp)) stop("Number of requested samples is higher than the total number of resampled entries. Check requested values again.")
      .tmp = .tmp[sample(1:nrow(.tmp), resample),]
    }

    # If the files we provided contained paths, we need to strip those out here
    strippedlogs = gsub(".*/","",logs)

    write.table(.tmp, file = paste0(out, "/", gsub(".Log.txt", "Log-NoHeader.txt", strippedlogs[l])), sep = "\t", col.names = T, row.names = F, quote = F)


    results[[gsub(".txt.Log.txt", "", logs[l])]] = .tmp
  }

  setwd(home)

  return(results)

}


#' @title Plot parameter values against index value.
#' @description Useful for visualizing parameter traces from MCMC analyses in
#' BayesTraits but can in principle be applied to any data.
#' @param files path or filename of the file a summary is required for.
#'     This can be a list of filenames - if so, all identical values in each
#'     table will be plotted on a single chart. This is useful for comparing estimates
#'     of the same parameter across multiple replicate,s for example.
#'     File must be a tab-delimited text file with column headers.
#'     This function is designed to be compatible with the output of [trimBTlog()].
#' @param cols A string or list of strings defining the column names or indices
#'     for which summary information is required. If unspecified, takes the
#'     value "all" which will return traces for the following columns
#'     (if present): Likelihood, Intercept, Slope(s)/Beta parameters,
#'     Variance, R-suqred, Number of local scalars (branch or node), and Lambda.
#' @param table if TRUE, then the function will accept a list of \code{data.frames}
#'     rather than a list of filenames.
#' @param logs if TRUE then the function will treat any filenames given as raw BayesTraits output (it will search for and remove any header info). Only compatible with table=FALSE and will not work with modified output files.
#' @param colours if provided, then should be a vector of colours for plotting - otherwise will be drawn randomly from rainbow.
#' @importFrom utils read.table
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines
#' @export
plotBTlog = function (files, cols = "all", table = T, logs = F, colours)
{
  if(logs)  cat("Reading output from raw log files...\n")
  # Check to see whether specified inputs are R objects or files
  if (table==T) {fil = files; files = names(files)} else {
    fil = list()
    for (f in 1:length(files))
    {
      if(!logs){
        fil[[f]] = utils::read.table(files[f], sep = "\t", header = T, stringsAsFactors = F)
      } else
        fil[[f]] = readBTlog(files[f])
      }
  }

  # If we have multiple files, create a legend.
  # Otherwise convert single file to list format.
  if(missing(colours))  colours = rainbow(n = length(files))
  if(!inherits(fil,"list")) {fil = list(fil)}  else{
    # Create a blank plot and add a colour legend IF more than one file
    plot(1, axes = F, main = "", bty = "n", type = "n",
         xlab = "", ylab = "")
    legend("center", legend = files, fill = colours)
    if(length(files) > 10) warning("Colours str specified randomly using rainbow() and may not be clearly distinct for more than 10 replicates.")
  }

  # Specify "iterations" as row indices
  its = seq(from = 1, to = max(unlist(lapply(fil, nrow))), by = 1)
  if (length(cols) == 1)
    if (cols == "all")
      cols = colnames(fil[[1]])[grepl("Lh|^Alpha|^Beta|Sigma|Var$|R.2|Local|Lambda",colnames(fil[[1]]))]


  # Loop through each of the specified columns of interest
  for (col in cols) {

    # Identify plot bounds using the full range of observed values
    bound = c(min(unlist(lapply(fil, function(x)min(x[,col])))),
              max(unlist(lapply(fil, function(x)max(x[,col])))))

    # Create a blank plot
    plot(1, xlim = c(min(its), max(its)), ylim = bound,
         type = "n", las = 1, xlab = "Iteration", ylab = col,
         bty = "l", main = col)

    # Add the values to the plot
    lapply(1:length(fil), function(x) lines(fil[[x]][, col] ~
              its[1:nrow(fil[[x]])], col = colours[x]))
  }
}

#' @title Extract variance from BayesTraits output.
#' @description This will work with raw output.
#' @param lf The path and filename of the log file the variance is wanted for.
#' @return List with elements:
#' \item{BG}{a vector of the extracted phylogenetic variance}
#' \item{dBG}{computed density distribution}
#' \item{dLogBG}{computed logged density distribution}
#' @importFrom stats density
#' @export
# Get variance from a BT log file
varBTlog = function(lf){
  skip = (which(grepl("\tTree No", readLines(lf))))-1
  log = read.table(lf, sep = "\t", header = T, stringsAsFactors = F, comment.char = "*", skip = skip, fill = T)
  BG = log[,which(grepl("Sigma|^Var$", colnames(log)))]
  return(list(BG = BG, dBG = density(BG), dLogBG = density(log10(BG))))
}



#' @title Function to summarize BayesTraits log files
#' @description This function will produce average parameter values (mean, median, mode) as well as other summary statistics.
#'     It will return averages, effective sample sizes, pMCMC values, and ranges.
#'     It ideally works with the output from [trimBTlog()].
#'     Can work with all sorts of outputs, but will need to manually specify column names for non Bayes Traits tables.
#' @param file data.frame or file path of the BayesTraits output.
#'    This must be one of either (see also the 'table' argument):
#'      - a file path for a tab-delimited text file specifying the output columns from a BayesTraits run.
#'      - a name of an  R data-frame.
#' @param cols A string or list of strings defining the columns for which summary information is required.
#'     If unspecified, takes the value "all" which summarizes a small pre-defined list of columns including:
#'     Likelihood, Alpha, Beta estimates, Variance, R-squared, Local transforms, and Lambda.
#' @param tradeoffs Boolean operator that if true, plots trade-offs between all parameters specified in cols.
#' @param table Boolean operator that if true, accepts an R data.frame object as input. If False, the file argument is interpreted as a character-string file path to the trimmed log file output from BayesTraits.
#' @param input Optional parameter. Can take two values. If simply specified as TRUE, this assumes that 'file' is a filename (i.e. table = F) and will search for a file with the name as follows: gsub(".Log.txt", "", file). Otherwise, can take a specified input file name. In either case, when specified, the function will modify the column names of the output table to reflect the input data file (e.g. Beta.1 will become "BodyMass"). Note that if the order of columns of the output or input files have been modified in any way this will not produce desirable output.
#' @param name The name which the output will be saved under.
#' @param logs Boolean operator that if true, means that the output file is in raw form as output directly from BayesTraits. Otherwise, it assumes the log file has been stripped of the header rows and is a simple tab-delimited file.
#' @importFrom grDevices pdf dev.off
#' @importFrom coda effectiveSize
#' @importFrom stats median density quantile
#' @return R data.frame object with the summary information:
#' \item{Parameter}{Name of the parameter of interest}
#' \item{Minimum}{Minimum parameter value}
#' \item{Mean}{Mean parameter value}
#' \item{Median}{Median parameter value}
#' \item{Mode}{Mode parameter value}
#' \item{Maximum}{Maximum parameter value}
#' \item{PropLess0}{The proportion of the posterior distribution of the parameter that is less than zero.}
#' \item{pmcmc}{The proportion of the posterior distribution of the parameter that crosses zero i.e. in either direction.}
#' @export


summarizeBTlog = function (file, cols = "all", tradeoffs = T, table = T, name = "Summary", input=T, logs = T) {
  out = NULL
  if (table) {
    fi = file
    file = name
  }  else
    if(!logs) {fi = read.table(file, sep = "\t", header = T)} else {fi = readBTlog(file)}


  if (length(cols) == 1)
    if (cols == "all"){
      cols = colnames(fi)[grepl("Lh|^Alpha|^Beta|Sigma|Var$|R.2|Local|Lambda",colnames(fi))]

      reorder=TRUE
    } else reorder = FALSE

    for (col in cols) {
      .tmp = fi[, col]
      l0 = length(.tmp[.tmp < 0])/length(.tmp)
      pmcmc = ifelse(l0 > 0.5, 1 - l0, l0)
      mod = density(.tmp)$x[which.max(density(.tmp)$y)]
      out = rbind(out, data.frame(Parameter = col, Minimum = min(.tmp),
                                  Mean = mean(.tmp), Median = median(.tmp), Mode = mod,
                                  Maximum = max(.tmp), PropLess0 = l0, pMCMC = pmcmc, q05 = quantile(.tmp, 0.05), q95 = quantile(.tmp, 0.95),
                                  stringsAsFactors = F))
    }
    ess = effectiveSize(fi[, cols])
    ess = data.frame(Parameter = names(ess), ess = ess)
    out = merge(out, ess)
    if (tradeoffs) {
      plot(fi[, cols])
    }

    out$Parameter[grepl("Beta", out$Parameter)]

    nums = as.numeric(as.character(gsub("Beta.", "", out$Parameter[grepl("Beta", out$Parameter)])))
    names(nums) = out$Parameter[grepl("Beta", out$Parameter)]
    nums = sort(nums)

    rownames(out) = out$Parameter
    out = rbind(out[!grepl("Beta", out$Parameter),], out[names(nums),])

    if(!missing(input)){
      if(input == T) input = gsub(".Log.txt", "", file) else input = input
      cat("Taking column headers from input file. Note this will not work if column order has been modified in any way.\n\n")
      inf=read.table(input, sep = "\t", header = T, stringsAsFactors = F)
      out$Parameter[grepl("Beta", out$Parameter)] = colnames(inf)[3:ncol(inf)]
      out$Parameter[grepl("Alpha", out$Parameter)] = colnames(inf[2])
    }

    # Restructure the table
    if(reorder){
      cols = c("Alpha", rownames(out)[grepl("Beta", rownames(out))], "Lh", "Var", "R.2")
      cols = c(cols, rownames(out)[!rownames(out) %in% cols])
      rows = c("Parameter", "Median", "pMCMC")
      rows = c(rows, colnames(out)[!colnames(out) %in% rows])
      out = t(out[cols,rows])

    }


    write.table(out, file = paste0(name, ".txt"), sep = "\t", col.names = F, row.names = T, quote = F)
    return(out)}

#' @title Calculate predictions from BayesTraits regressions.
#' @description Function to calculate predictions from BayesTraits input/output files. Incompatible with other programs.
#' @param input A character vector defining the original input file (must have column headers).
#' Note that this can be a modified version of the original input e.g. if one wants to fix a value at a mean - as long as columns are retained in identical order and format.
#' @param output an optional argument. If defined, a character vector locating the original output file (unmodified).
#' If unspecified, the function will search for the original input file as if it had been run through BayesTraits, appending .Log.txt.
#' @param plot Optional argument. If TRUE, the function will generate a prediction plot. For complicated datasets this might produce undesirable results, and so the behaviour defaults to FALSE.
#' This should really only be used in the case of a simple linear regression (i.e. a single continuous predictor and a single continuous response).
#' @return A matrix of predicted values for terminal taxa (if original input was MCMC, one iteration per column)
#' @export
predBTlog = function(input, output = NULL, plot = F){

  # Specify output (if not specified)
  if(is.null(output)) output = paste0(input, ".Log.txt")

  # Read in tables
  input = read.table(input, sep = "\t", header = T, stringsAsFactors = F)
  output = read.table(file = output, sep = "\t", header = T, skip = which(grepl("Sites:", readLines(output))), stringsAsFactors = F)

  # Extract Xs and Bs
  Xs = input[,3:ncol(input), drop = F]
  Bs = output[,grepl("^Alpha|^Beta", colnames(output)), drop = F]

  # Predict across all Xs
  preds = apply(Bs, 1, function(y)unlist(apply(Xs, 1, function(x)sum(x*y[2:length(y)]) + y[1])))

  if(plot == T){
    plot(input[,2] ~ input[,3], type = "n")
    apply(preds,2,function(x)lines(x~input[,3], col = makeTransparent("black", alpha = 0.1)))}

  rownames(preds) = input[,1]
  return(preds)

}

#' @title Extract marginal likelihood from a stones file
#' @description Function to extract the marginal likelihood from the raw stepping-stone sampling output.
#' @param input A character vector defining the stones file exactly as output by BayesTraits.
#' @return A numeric value defining the log marginal likelihood of the output model.
#' @export
extractML = function(input){

  # read in raw file
  raw = readLines(input)

  # Get the last line
  Lhline = raw[length(raw)]

  # Get the value
  Lh = as.numeric(gsub(".+\t", "", Lhline))

  return(Lh)

}
