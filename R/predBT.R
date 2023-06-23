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
predBT = function(input, output = NULL, plot = F){

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
