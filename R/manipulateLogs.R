#' @title Read in BayesTraits Log Files (raw)
#' @description Removes the header from BayesTraits output files.
#'     Can only work with log files (.txt.Log.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not work with output from any other program or model.
#' @param file Takes a single string defining the direct path to a single log file.
#' @return R data.frame object that contain the trimmed log file.
#' @export

readlog = function(file){
  skip = (which(grepl("\tTree No", readLines(file)))) - 1
  log = read.table(file, sep = "\t", header = T, stringsAsFactors = F,
                   comment.char = "*", skip = skip, fill = T)
  return(log)
}


#' @title Trim header from BayesTraits log files.
#' @description Removes the  header from BayesTraits output files.
#'     Can only work with log files (.txt.Log.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not work with output from any other program or model.
#' @param dir Directory of files for which we want trimmed output. Defaults to the current working directory.
#' @param out The name of the folder to which we want output saved to.
#' @param pat If no files are specified, the function will use the string defined by pat to search for files in dir.
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
#'     If this is defined, dir and pat are ignored.
#' @return R data.frame objects that contain the trimmed log files.
#' @export


trimmedlogs = function(dir = ".", out = "TrimmedLogFiles", pat = ".txt.Log.txt", start = "\tTree No", tail = NA, resample = NA, files= NULL){

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
    .tmp = readlog(logs[l])

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
#'     This function is designed to be compatible with the output of [trimmedlogs()].
#' @param cols A string or list of strings defining the column names or indices
#'     for which summary information is required. If unspecified, takes the
#'     value "all" which will return traces for the following columns
#'     (if present): Likelihood, Intercept, Slope(s)/Beta parameters,
#'     Variance, R-suqred, Number of local scalars (branch or node), and Lambda.
#' @param table if TRUE, then the function will accept a list of \code{data.frames}
#'     rather than a list of filenames.
#' @param logs if TRUE then the function will treat any filenames given as raw BayesTraits output (it will search for and remove any header info). Only compatible with table=FALSE and will not work with modified output files.
#' @importFrom utils read.table
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines
#' @export
plotTraces = function (files, cols = "all", table = T, logs = F)
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
        fil[[f]] = readlog(files[f])
      }
  }

  # Specify "iterations" as row indices
  its = seq(from = 1, to = max(unlist(lapply(fil, nrow))), by = 1)
  if (length(cols) == 1)
    if (cols == "all")
      cols = colnames(fil[[1]])[grepl("Lh|^Alpha|^Beta|Sigma|Var$|R.2|Local|Lambda",
                                      colnames(fil[[1]]))]

  # Create a blank plot and add a colour legend
  colours = rainbow(n = length(files))
  plot(1, axes = F, main = "", bty = "n", type = "n",
       xlab = "", ylab = "")
  legend("center", legend = files, fill = colours)
  if(length(files) > 10) warning("Colours str specified randomly using rainbow() and may
                                 not be clearly distinct for more than 10 replicates.")

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
getVar = function(lf){
  skip = (which(grepl("\tTree No", readLines(lf))))-1
  log = read.table(lf, sep = "\t", header = T, stringsAsFactors = F, comment.char = "*", skip = skip, fill = T)
  BG = log[,which(grepl("Sigma|^Var$", colnames(log)))]
  return(list(BG = BG, dBG = density(BG), dLogBG = density(log10(BG))))
}



#' @title Function to summarize BayesTraits log files
#' @description This function will produce average parameter values (mean, median, mode) as well as other summary statistics.
#'     It will return averages, effective sample sizes, pMCMC values, and ranges.
#'     It ideally works with the output from trimmedlogs().
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
#' @param name The name which the output will be saved under.
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

summarizeBT = function (file, cols = "all", tradeoffs = T, table = T, name = "Summary") {
  out = NULL
  if (table) {
    fi = file
    file = name
  }
  else fi = read.table(file, sep = "\t", header = T)
  if (length(cols) == 1)
    if (cols == "all")
      cols = colnames(fi)[grepl("Lh|^Alpha|^Beta|Var$|R.2|Local|Lambda",
                                colnames(fi))]
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
    pdf(paste0(name, "_Tradeoffs.pdf"), useDingbats = F,
        height = 20, width = 20)
    plot(fi[, cols])
    dev.off()
  }

  out$Parameter[grepl("Beta", out$Parameter)]

  nums = as.numeric(as.character(gsub("Beta.", "", out$Parameter[grepl("Beta", out$Parameter)])))
  names(nums) = out$Parameter[grepl("Beta", out$Parameter)]
  nums = sort(nums)

  rownames(out) = out$Parameter
  out = rbind(out[!grepl("Beta", out$Parameter),], out[names(nums),])

  write.table(out, file = paste0(name, ".txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  return(out)}

