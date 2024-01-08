#' @title Format VR block from BayesTraits log file
#' @description Helper function used within readVR to format the VR block of a
#' variable rates output file into an R-friendly data.frame.
#' @param VRblock Text string including all rows of the VR output relevant to the
#' specification and definition of rate scalars. Created within the \code{readVR}
#' function.
#' @param colstart Numeric value specifying the column number where scalar information
#' begins to be recorded in the VRblock output. This varies between different types
#' of VR model and is internally identified within the \code{readVR} function.
#' @returns A \code{data.frame} recording detailed information about the rate
#' scalars acting at each iteration of the variable rates model.
#' @export
#'
tidyVRblock <- function(VRblock,colstart){

  # Print output
  cat("\tConverting output to tabular format.\n")

  # Extract column names
  .vrlognames = VRblock[[1]]
  .vrlogblock = VRblock[-1]

  # Tidy up and split names of VR log block
  .vrlognames_1 = gsub(" ", "",.vrlognames[1:(colstart-1)])
  .vrlognames_2 = gsub(" |/", "",.vrlognames[colstart:length(.vrlognames)])

  # Info for progress updates
  itnumber=seq(from = 1, to = length(.vrlogblock), by = length(.vrlogblock)/10)
  names(itnumber) = seq(from = 0, to = 90, by = 10)

  # Loop through vr log block and add names
  vrlogtable = NULL
  for(i in 1: length(.vrlogblock)){
    # Print some output
    if(i %in% itnumber)
      cat("\t\t - ", names(itnumber)[itnumber == i], "%\n")

    if((length(.vrlogblock[[i]])-(colstart-1)) %% length(.vrlognames_2) != 0)
      stop("vr logfile column headers not as expected.
                   Function currently only compatible with variable rates.")
    .names = rep(.vrlognames_2,(length(.vrlogblock[[i]])-(colstart-1)) / length(.vrlognames_2))
    names(.vrlogblock[[i]]) = c(.vrlognames_1, .names)
    # }
    # Tidy up
    .tmp = .vrlogblock[[i]]
    .firsthalf = data.frame(as.list(.tmp[1:(colstart-1)]))
    if(length(.tmp) == length(.vrlognames_1)) {
      .secondhalf = data.frame(CreatIt= NA, NodeBranch=NA,PartID=NA,Scaler=NA)
      } else  .secondhalf = data.frame(split(.tmp[colstart:length(.tmp)], names(.tmp[colstart:length(.tmp)])))
    .vrdf = cbind(.firsthalf,.secondhalf)
    vrlogtable = rbind(vrlogtable,.vrdf)

  }

  # Return
  cat("\t\t - done")
  return(vrlogtable)
}

#' @title Read in BayesTraits VR Log Files (raw)
#' @description Reads the output file from a variable rates analysis in BayesTraits
#'     for further processing in an R workspace. Will only work with VR log files
#'     (.VarRates.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not (currently) work with output from any other program or model.
#' @param vrfile Takes a single string defining the direct path to a single log file.
#' @returns A list with two elements:
#' \itemize{
#' \item brInfo - A \code{data.frame} giving information about each branch/node
#' in the analysis. Each branch is identified by the numeric value in "branch". For
#' samples of trees, this value is equivalent to the "partition ID" used internally
#' by BayesTraits. "NoDes" and "Des" record the number and names of all descendant
#' taxa from each branch. The "freq" and "probability" columns tell us the frequency
#' and probability that this branch occur in the sample are only relevant for
#' analyses over samples of trees. The "length" records the branch length and
#' is left blank for tree samples.
#' \item VRlog - A data.frame tabulating the VR output from BayesTraits and
#' containing the following columns: "It", "Lh", "Lh.Prior", "NoPram", "Alpha",
#' "Sigma.2", "AlphaScalePrior", "CreateIt", "NodeBranch", "NodeID", "Scalar". If
#' the original analysis was done on a sample of trees, "NodeID" is replaced by
#' "Partition ID" and the tree number used at each iteration is recorded.
#' }
#' @export
#' @importFrom utils head
readVR <- function(vrfile){

  # Read full input file
  vrraw = readLines(vrfile)

  # Identify vr log block
  .vrlogstart = which(grepl("^It", vrraw))
  .vrlog = strsplit(vrraw[(.vrlogstart):length(vrraw)], split = "\t")

  # Print output
  cat("Extracting VR information from file...\n")

  # If we have a sample of trees, the structure of the output is different
  if(grepl("Part", vrraw[1])){

    # Flag up that this is multi-topology
    ntrees = Inf

    # Number of node/branch combinations
    nbranches = as.numeric(unlist(strsplit(vrraw[1],split = "\t"))[2])

    # Extract branch information block
    brInfoblock = vrraw[2:nbranches+2]
    brInfolist = sapply(brInfoblock, function(x)strsplit(x,split = "\t"), USE.NAMES=F)
    brInfotable = data.frame(t(sapply(brInfolist, utils::head, 5)), row.names = NULL)
    colnames(brInfotable) = c("branch", "freq", "probability", "length", "NoDes")

    # Tidy up descendants
    brDes = sapply(brInfolist, "[", -c(1:5))
    names(brDes) = brInfotable$branch

    # Identify where the columns of scalar information begin
    .colstart = 9

  }else {

    # This is only a single topology
    ntrees =1

    # Number of taxa
    ntax = as.numeric(vrraw[1])

    # Identify taxa definitions
    txdef = read.table(vrfile, skip = 1, nrows = ntax,
                       sep = "\t",  header = F, col.names = c("TaxNo", "Tax"))

    # Number of nodes
    nnode = as.numeric(vrraw[ntax+2])

    # Identify node definition block and extract info
    nddef = vrraw[(ntax+3):(nnode+ntax+2)]
    nddef = lapply(nddef, function(x)unlist(strsplit(x, split = "\t")))
    branchIDs = sapply(nddef,"[[", 1)

    # Pull out branch info table
    brInfotable = data.frame(t(sapply(nddef, head, 3)), row.names = NULL)
    colnames(brInfotable) = c("branch", "length", "NoDes")
    brInfotable[,c("freq", "probability")] = NA
    brInfotable = brInfotable[,c("branch", "length", "freq", "probability", "NoDes")]

    # Pull out branch descendants
    brDes = sapply(nddef, "[", -c(1,2,3))
    brDes = lapply(brDes,function(y)unlist(lapply(y,function(x)txdef$Tax[match(x,txdef$TaxNo)])))
    names(brDes) = branchIDs


    # Identify where the columns of scalar information begin
    .colstart = 8

  }


  # Tidy VR block
  vrlogtable = tidyVRblock(.vrlog, .colstart)

  # Add descendant information to table
  brInfotable$Des = unlist(sapply(brDes,function(x)paste0(sort(x), collapse = ",")))

  return(list(brInfo = brInfotable, VRlog = vrlogtable, ntrees = ntrees))
}


