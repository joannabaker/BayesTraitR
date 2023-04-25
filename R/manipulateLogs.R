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

