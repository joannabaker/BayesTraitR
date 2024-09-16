## TODO: Clarify argument descriptions for all functions.
## TODO: Check compatibility with outputs with equal trees.
## TODO: link trees is way too slow with big trees, fix please!


#' @title Read in BayesTraits VR Log Files (raw)
#' @description Reads the output file from a variable rates analysis in BayesTraits
#'     for further processing in an R workspace. Will only work with VR log files
#'     (.VarRates.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#' @param vrfile Takes a single string defining the direct path to a single variable rates output file.
#' @param burnin A value specifying the number of rows to remove from the beginning of the sampled chain. Default value is 0 and so the entire chain will be read in.
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
readVR = function(vrfile){

  # Read full input file
  vrraw = readLines(vrfile)

  # Identify vr log block
  .vrlogstart = which(grepl("^It", vrraw))
  .vrlog = strsplit(vrraw[(.vrlogstart):length(vrraw)], split = "\t")

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

  # Remove burnin
  .vrlogblock = .vrlogblock[(burnin+1):length(.vrlogblock)]

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

  # If multiple trees, change the column names for consistency
  if(ntrees !=1) {
    VRblock[,"Node ID"] = VRblock[,"Part ID"]
    VRblock[,"Part ID"] = NULL
  }

  return(list(brInfo = brInfotable, VRlog = VRblock, ntrees = ntrees))
}


#' @title Link branches between processed VR output and an input tree
#' @description Identify corresponding branches between variable rates output and a user-defined input tree.
#' @param VRout An R object defining the output of the read VR function.
#' @param tree A phylogenetic tree object of class \code{phylo}. The output will link rates to the given tree. That is, every rate scalar is defined on the basis of a branch. If that branch does not exist in the defined tree, the rate scalar will instead be applied to the branch which terminates in the most recent common ancestor of all taxa that descended from the original defined branch. This is designed to summarize the results of a multi-topology variable rates output. The input tree must contain all taxa in the original tree used to run the variable rates model.
#' @keywords internal
#' @noRd
#' @importFrom ape read.nexus
#' @importFrom phytools findMRCA
#' @importFrom tidyr unnest

linkbranches = function(VRout, tree){

  # Split VRout into branch info and branch output
  VRbrInfo = VRout$brInfo
  VRlog = VRout$VRlog

  # Tabulate the tree
  treetab = tabulatetree(tree)

  # Modify column names to fit with the output of readVR
  colnames(treetab)[c(5, 7)] = c("Des", "DesBranches")

  # Add a row for the root
  treetab = rbind(c((length(tree$tip.label)+1),"NA", 0, -1, paste0(sort(tree$tip.label),collapse=","), Ntip(tree),paste0(0:nrow(treetab), collapse = ",")), treetab)

  # Get descendant lists for VRinfo
  deslist = lapply(VRbrInfo$Des,function(x)unlist(strsplit(x,split=",")))
  names(deslist) = VRbrInfo$branch

  # Remove root from list
  if(names(deslist)[1] == "0") deslist = deslist[2:length(deslist)]

  # Split tips and nodes
  des_tips = deslist[which(unlist(lapply(deslist,length)) < 2)]
  des_nodes = deslist[which(unlist(lapply(deslist,length)) > 1)]

  # Check - need trees with the same taxa in
  if(any(lapply(des_nodes, function(x)all(x %in% tree$tip.label)) == FALSE))
    stop("All taxa in the rates tree must be found in the input tree.")

  # identify the MRCA of each of the nodes in the new tree
  mrcas = sapply(des_nodes, function(x)treetab$Rbranch[treetab$DescNode == findMRCA(tree, tips = x)])

  # identify the tip label of each of the tips in the new tree
  tips = sapply(des_tips, function(x)treetab$Rbranch[treetab$Des == x])

  # Combine them back into a reference list
  reference = c(mrcas, tips)

  # Modify all branch references in the log
  VRlog[,"Rbranch"] = reference[match(VRlog[,"Node ID"], names(reference))]

  # Output the results
  VRout = list(brInfo = VRbrInfo, VRlog = VRlog, treetab = treetab, tree = tree)
  return(VRout)

}


#' @title Process BayesTraits VR log file
#' @description Process the output files from a variable rates analysis in BayesTraits as read in by readVR.
#' @param VRout An R object defining the output of the read VR function with branches linked to R.
#' Designed for use with the output of linkbranches. This will fail if used without running linkbranches first.
#' @returns A list with three elements:
#' \itemize{
#' \item summary - A \code{data.frame} giving information about each branch/node
#' in the analysis. Each branch is identified by the numeric value in "branch". For
#' samples of trees, this value is equivalent to the "partition ID" used internally
#' by BayesTraits. "NoDes" and "Des" record the number and names of all descendant
#' taxa from each branch. The "freq" and "probability" columns tell us the frequency
#' and probability that this branch occur in the sample are only relevant for
#' analyses over samples of trees. The "length" records the branch length and
#' is left blank for tree samples. The "DesBranches" tells us all branches that descend from
#' the branch of interest. The median and mean scalar are the median and mean rate
#' scalars acting on that branch across the entire posterior distribution. n_scaled
#' tells us how many iterations that branch had a scalar acting on it, and n_origin tells
#' us how many iterations that branch had a scalar ORIGINATING on it (either as the base of
#' a node scalar or as its own branch scalar).
#' \item posterior - A \code{data.frame} providing the posterior distribution of rate scalars acting
#' on every branch (column) at each iteration (row).
#' \item tree - The tree will be returned for use in tree-scaling and visualization etc.
#' }
#' @keywords internal
#' @noRd
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyr unnest
processVR = function(VRout){

  # First identify iterations and branches
  uiterations = unique(VRout$VRlog$It)
  ubranches = VRout$treetab$Rbranch

  # Create a table for the posteriors
  posterior <- matrix(1, nrow = length(uiterations), ncol = length(ubranches),
                      dimnames = list(uiterations, ubranches))

  # OK, now we want to split nodes and branches
  nodes <- VRout$VRlog %>% dplyr::filter(NodeBranch == "Node") %>% select(It, Scaler, Rbranch)
  nodes$Rbranch = sapply(nodes$Rbranch, function(x)unlist(strsplit(VRout$treetab$DesBranches[VRout$treetab$Rbranch == x], split = ",")))
  #  nodes$Rbranch = sapply(nodes$Rbranch, function(x)unlist(strsplit(VRout$treetab$DesBranches[VRout$treetab$Rbranch == x], split = ","))[-1])


  # Expand the nodes so that we have a row per branch
  nodes_expanded = nodes %>%
    unnest(Rbranch)
  nodes_expanded = as.data.frame(nodes_expanded)

  # Now format the branch table in the same way
  branches = VRout$VRlog %>% dplyr::filter(NodeBranch == "Branch") %>% select(It, Scaler, Rbranch)

  # Combine the two
  allscalars = rbind(nodes_expanded,branches)

  # Now aggregate the scalar values per iteration (no duplicates)
  allscalars$Scaler=as.numeric(allscalars$Scaler)
  aggregated_raw <- allscalars %>%
    dplyr::group_by(It,Rbranch) %>%
    dplyr::summarise(Scaler=prod(Scaler),.groups='drop')

  # Set a chunk size for processing to avoid memory issues
  chunk_size = 5000

  # Calculate the number of chunks needed
  num_chunks = ceiling(nrow(aggregated_raw)/ chunk_size)

  # Initialize posterior matrix index
  posterior_index = matrix(NA, nrow = nrow(aggregated_raw), ncol = 2)
  posterior_index[, 1] = match(aggregated_raw$It, uiterations)
  posterior_index[, 2] = match(aggregated_raw$Rbranch, ubranches)

  # Update the posterior matrix chunk by chunk
  # Reduce chunk size if R runs out of memory
  # Increase chunk size to speed up operation
  for (chunk in 1:num_chunks) {
    start_idx = (chunk - 1) * chunk_size + 1
    end_idx = min(chunk * chunk_size, nrow(aggregated_raw))

    # Subset of posterior index
    index_subset <- posterior_index[start_idx:end_idx, ]

    # Subset of aggregated_raw
    aggregated_subset <- aggregated_raw[start_idx:end_idx, ]

    # Update posterior matrix
    posterior[index_subset] <- aggregated_subset$Scaler
  }

  # Add summary info
  summarytab = VRout$treetab
  summarytab$medianscalar = unlist(apply(posterior, 2, median))
  summarytab$meanscalar = unlist(apply(posterior, 2, mean))
  summarytab$n_scaled = apply(posterior,2,function(x)(length(x[x>1])+length(x[x<1])))

  # Return
  return(list(ratesummary = summarytab, posterior = posterior, tree = VRout$tree))
}



#' @title Summarize BayesTraits VR Log Files (from raw output)
#' @description Reads the output file from a variable rates analysis in BayesTraits
#'     and produces summarized output. Will only work with VR log files
#'     (.VarRates.txt) as directly output from BayesTraits models 1-4, 7 and 9.
#'     Will not (currently) work with output from any other program or model.
#' @param vrfile Takes a single string defining the direct path to a single variable rates log file.
#' @param tree A phylogenetic tree object of class \code{phylo}. The output will link rates to the given tree. That is, every rate scalar is defined on the basis of a branch. If that branch does not exist in the defined tree, the rate scalar will instead be applied to the branch which terminates in the most recent common ancestor of all taxa that descended from the original defined branch. This is designed to summarize the results of a multi-topology variable rates output. The input tree must contain all taxa in the original tree used to run the variable rates model.
#' @param forcedequaltrees Boolean operator, default = \code{FALSE}. If \code{forcedequaltrees = true}, then the model will be run iteratively over a pre-defined sample of trees.In this situation, we will have duplicate iterations (one for each tree) and thus multiple posteriors. We need to let the post-processor know if this is the case.
#' @param burnin A value specifying the number of rows to remove from the beginning of each log file. Default value is 0 and so the entire chain will be used to calculate rates.
#' @returns A list with three elements:
#' \itemize{
#' \item summary - A \code{data.frame} giving information about each branch/node
#' in the analysis. Each branch is identified by the numeric value in "branch". For
#' samples of trees, this value is equivalent to the "partition ID" used internally
#' by BayesTraits. "NoDes" and "Des" record the number and names of all descendant
#' taxa from each branch. The "freq" and "probability" columns tell us the frequency
#' and probability that this branch occur in the sample are only relevant for
#' analyses over samples of trees. The "length" records the branch length and
#' is left blank for tree samples. The "DesBranches" tells us all branches that descend from
#' the branch of interest. The median and mean scalar are the median and mean rate
#' scalars acting on that branch across the entire posterior distribution. n_scaled
#' tells us how many iterations that branch had a scalar acting on it, and n_origin tells
#' us how many iterations that branch had a scalar ORIGINATING on it (either as the base of
#' a node scalar or as its own branch scalar).
#' \item posterior - A \code{data.frame} providing the posterior distribution of rate scalars acting
#' on every branch (column) at each iteration (row).
#' \item tree - The input tree will be returned for use in tree-scaling and visualization etc.
#' }
#' @examples
#' # Here we define a file based on one of the examples included with this package
#' VRlogfilepath <- system.file("extdata", "MammalBody_VR-001.txt.VarRates.txt")
#' summarizeVR(vrfile = VRlogfilepath, tree = Mammal_trees)
#' @export
summarizeVR <- function(vrfile, tree, forcedequaltrees = F, burnin = 0){
  # Read log
  cat("Extracting VR information from file...\n")
  VRout = readVR(vrfile, burnin = burnin)

  # Link to tree
  cat("\tLinking to tree...\n")
  VRout = linkbranches(VRout, tree)

  if(forcedequaltrees){
    VRout$VRlog$It = paste0(VRout$VRlog$It, "_", VRout$VRlog$`Tree No`)
  }

  # Link to tree
  cat("\t\tSummarizing Rates...\n")
  summarized = processVR(VRout)
  return(summarized)
}


#' @title Create a stretched tree from variable rates output
#' @description Creates a stretched tree (mean, median, or full sample) from the output of the variable rates model.
#' @param VRsummary The name of the object created by [summarizeVR]
#' @param frequency A percentage ranging between 0-100 specifying which branches to scale. Branches will only be scaled if they have been scaled in a greater percentage of the posterior than is specified by this parameter.This parameter is ignored if type = "sample".
#' @param magnitude A value specifying wich branches to scale. Branches will only be scaled if they are above the value specified by magnitude (or below 1/magnitude for rate decreases). This parameter is ignored if type = "sample".
#' @param type A character string defining what type of scaled tree to return. Can be one of "median", "mean", or "sample". If sample, the full posterior of scaled trees will be returned.
#' @returns A tree where branch lengths are stretched by the mean or median scalar - or a list of trees stretched by the rate scalars.
#' @export
scaleTree = function(VRsummary, frequency = 0, magnitude = 1, type = "median"){

  # Check type of scaled tree desired
  if(!type %in% c("median", "mean", "sample")) stop("type argument must be one of median, mean, or sample.")

  # Get the tree
  tree = VRsummary$tree

  # If sample, create sample
  if(type == "sample"){
    trees = list(time = tree)
    for(i in 2:nrow(VRsummary)){
      tree1 = tree
      tree1$edge.length = tree1$edge.length * VRsummary$posterior[i,2:ncol(VRsummary$posterior)]
      trees[[i]] = tree1
    }
    class(trees) = "multiPhylo"
    return(trees)
  }

  # get the scalar
  scalar = VRsummary$ratesummary[,paste0(type, "scalar")]

  # calculate the frequency
  freq = VRsummary$ratesummary$n_scaled/nrow(VRsummary$posterior)*100

  # If less than the specified frequency, set scalar to 1
  scalar[freq < frequency] = 1

  # Identify the absolute magnitude
  scalar[scalar < magnitude & scalar > 1] = 1
  scalar[scalar < 1 & scalar < 1/magnitude]  = 1

  # Scale the tree
  tree$edge.length = tree$edge.length * scalar[-1]

  return(tree)
}

