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
#' @keywords internal
#' @noRd
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

  }

  # Tidy VR block
  .vrlognames = .vrlog[[1]] # Identify column headers
  .vrlogblock = .vrlog[-1] # Remove them from the block

  # Split the log block into two parts:
  # First - per-iteration details
  .vrheader = lapply(.vrlogblock, function(x)x[1:(length(.vrlognames)-4)])

  # Second - scalar details
  .vrlogblock = lapply(.vrlogblock, function(x)x[(length(.vrlognames)-3):length(x)])

  # Convert scalars to a matrix
  .vrlogblock = lapply(.vrlogblock,matrix,ncol=4,byrow=T)

  # Add the per-iteration details to the scalar matrix
  .vrlogblock = lapply(1:length(.vrlogblock),function(x)cbind(matrix(rep(.vrheader[[x]], each = nrow(.vrlogblock[[x]])), ncol = length(.vrheader[[x]]), byrow = F),.vrlogblock[[x]]))

  # Bind them together
  VRblock = do.call(rbind,.vrlogblock)

  # Convert to a data frame
  VRblock = as.data.frame(VRblock)
  colnames(VRblock) = .vrlognames

  # Adjust column name
  colnames(VRblock)[ncol(VRblock)] = "NodeBranch"

  # Add descendant information to table
  brInfotable$Des = unlist(sapply(brDes,function(x)paste0(sort(x), collapse = ",")))

  # Now also add all branches descending to the table, too
  alldes = lapply(brInfotable$Des, function(x)unlist(strsplit(x,split=",")))
  brInfotable$DesBranches = unlist(lapply(alldes,function(y)paste0(brInfotable$branch[which(unlist(lapply(alldes,function(x)all(x %in% y))))], collapse = ",")))

  # If multiple trees, change the column names for consistency
  if(ntrees !=1) {
    VRblock$NodeID = VRblock[,"Part ID"]
    VRblock[,"Part ID"] = NULL
  }


  return(list(brInfo = brInfotable, VRlog = VRblock, ntrees = ntrees))
}
#'
#'
#'
#' #' @title Process BayesTraits VR log file
#' #' @description Process the output files from a variable rates analysis in BayesTraits
#' #'     as read in by readVR.
#' #' @param VRout An R object defining the output of the read VR function.
#' #' @param treefile An optional parameter giving the name of a treefile as a string.
#' #' If defined, the output will link rates to the given tree. That is, every rate scalar
#' #' is defined on the basis of a branch. If that branch does not exist in the defined tree, the
#' #' rate scalar will instead be applied to the branch which terminates in the most recent
#' #' common ancestor of all taxa that descended from the original defined branch. This is
#' #' designed explicitly to summarize the results of a multi-topology variable rates output.
#' #' @keywords internal
#' #' @noRd
#' #' @importFrom ape read.nexus
#' #' @importFrom phytools findMRCA
#' processVR = function(VRout, treefile){
#'
#'   # Print output
#'   cat("\tSummarizing scalars\n")
#'
#'   # Identify the unique iterations
#'   iters = unique(VRout$VRlog$It)
#'
#'   # If tree is specified, then we need to input our own branch table into the summary
#'   if(!missing(treefile)){
#'     cat("\tTree-file identified: linking to MRCAs.\n\tWARNING: This might be slow.\n")
#'     tree = read.nexus(treefile)
#'     treetab = tabulatetree(tree)
#'     colnames(treetab)[c(3,6)] = c("branch", "DesBranches")
#'     treetab = rbind(c("NA", (length(tree$tip.label)+1), 0, -1, paste0(sort(tree$tip.label),collapse=","), paste0(0:nrow(treetab), collapse = ",")), treetab)
#'   } else treetab = VRout$brInfo
#'
#'   # Create a summary table
#'   ratestable_summary = data.frame(VRout$brInfo)
#'   ratestable_summary[,c("median_scalar", "mean_scalar")] = NA
#'   ratestable_summary[,c("n_scaled", "n_origin" )] = 0
#'
#'   # Create a table for the posteriors
#'   ratestable_post = as.data.frame(matrix(1, ncol=nrow(ratestable_summary), nrow=length(unique(iters)), dimnames=list(NULL, ratestable_summary$branch)))
#'   rownames(ratestable_post) = iters
#'
#'   # Initialize the progress bar
#'   progress_bar = txtProgressBar(min = 0, max = nrow(VRout$VRlog), style = 3)
#'
#'   # Now traverse the VRlog file and apply scalars to each branch at each iteration
#'   for(i in 1:nrow(VRout$VRlog)){
#'
#'     # Identify iteration
#'     iter = VRout$VRlog$It[i]
#'
#'     # If NO rate scalars identified, move on
#'     if(is.na(VRout$VRlog$Scaler[i])) next
#'
#'     # Identify branch
#'     if(!missing(treefile)){
#'
#'       nodeID = VRout$VRlog$NodeID[i]
#'       tax =  unlist(strsplit(VRout$brInfo$Des[VRout$brInfo$branch == nodeID], split = ","))
#'       newnode = findMRCA(tree, tax)
#'       branch = origin = treetab$branch[treetab$DescNode == newnode]
#'
#'
#'     }else {
#'       branch = origin = ratestable_summary$branch[which(ratestable_summary$branch == VRout$VRlog$NodeID[i])]
#'     }
#'
#'     # If this is a node scalar, this needs to be applied across all branches
#'     if(VRout$VRlog$NodeBranch[i] == "Node"){
#'       branch = unlist(strsplit(ratestable_summary$DesBranches[ratestable_summary$branch == branch], ","))
#'     }
#'
#'     # Apply the rate scalar
#'     ratestable_post[iter,branch] = ratestable_post[iter,branch]* as.numeric(VRout$VRlog$Scaler[i])
#'
#'     # Add n to the summary table
#'     ratestable_summary[branch, "n_scaled"] = ratestable_summary[branch, "n_scaled"] + 1
#'     ratestable_summary[origin, "n_origin"] = ratestable_summary[origin, "n_origin"] + 1
#'
#'     setTxtProgressBar(progress_bar, value = i)
#'
#'   }
#'
#'   # Calculate the mean and median rate scalars
#'   ratestable_summary$median_scalar = apply(ratestable_post,2, median)
#'   ratestable_summary$mean_scalar = apply(ratestable_post,2, mean)
#'
#'
#'   out = list(summary = ratestable_summary, posterior = ratestable_post)
#'   if(!missing(treefile)) out[["tree"]] = tree
#'
#'   # Close progress bar
#'   close(progress_bar)
#'
#'   return(out)
#'
#' }




#' #' @title Summarize BayesTraits VR Log Files (from raw output)
#' #' @description Reads the output file from a variable rates analysis in BayesTraits
#' #'     and produces summarized output. Will only work with VR log files
#' #'     (.VarRates.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#' #'     Will not (currently) work with output from any other program or model.
#' #' @param vrfile Takes a single string defining the direct path to a single log file.
#' #' @param treefile An optional parameter giving the name of a treefile as a string.
#' #' If defined, the output will link rates to the given tree. That is, every rate scalar
#' #' is defined on the basis of a branch. If that branch does not exist in the defined tree, the
#' #' rate scalar will instead be applied to the branch which terminates in the most recent
#' #' common ancestor of all taxa that descended from the original defined branch. This option is
#' #' designed explicitly to summarize the results of a multi-topology variable rates output.
#' #' @returns A list with two elements:
#' #' \itemize{
#' #' \item summary - A \code{data.frame} giving information about each branch/node
#' #' in the analysis. Each branch is identified by the numeric value in "branch". For
#' #' samples of trees, this value is equivalent to the "partition ID" used internally
#' #' by BayesTraits. "NoDes" and "Des" record the number and names of all descendant
#' #' taxa from each branch. The "freq" and "probability" columns tell us the frequency
#' #' and probability that this branch occur in the sample are only relevant for
#' #' analyses over samples of trees. The "length" records the branch length and
#' #' is left blank for tree samples. The "DesBranches" tells us all branches that descend from
#' #' the branch of interest. The median and mean scalar are the median and mean rate
#' #' scalars acting on that branch across the entire posterior distribution. n_scaled
#' #' tells us how many iterations that branch had a scalar acting on it, and n_origin tells
#' #' us how many iterations that branch had a scalar ORIGINATING on it (either as the base of
#' #' a node scalar or as its own branch scalar).
#' #' \item posterior - A \code{data.frame} providing the posterior distribution of rate scalars acting
#' #' on every branch (column) at each iteration (row).
#' #' }
#' #' @export
#' summarizeVR <- function(vrfile, treefile){
#'   VRout = readVR(vrfile)
#'   if(!missing(treefile)) VRsum = processVR(VRout, treefile) else VRsum = processVR(VRout)
#'   return(VRsum)
#' }
