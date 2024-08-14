#' @title Get all descendants of all nodes in a tree
#' @description Retrieves the descendants of every node (including terminals) in a tree.
#' @param tree a phylogenetic tree as an object of class \code{phylo}.
#' @return A \code{list} of vectors defining the taxa descending from every node of a tree.
#' @keywords internal
#' @export
branchdefs = function(tree){
  nodes = 1:(Ntip(tree)+Nnode(tree))
  des = lapply(nodes, function(x)paste0(sort(tipDescendants(tree, x)), collapse = ","))
  names(des) = nodes
  return(des)
}


#' @title Tabulate branch-level information for a tree.
#' @description Gets all information associated with every branch of a tree.
#' @param tree a phylogenetic tree as an object of class \code{phylo}.
#' @return A \code{data.frame} with \code{n} rows equal to the number of branches in the tree
#' and the following columns:
#'
#' * AncNode: the label for the node defining the start of the branch
#' * DescNode: the label for the node defining the end of the branch
#' * Rbranch: the branch ID
#' * length: the branch length
#' * des: The list of taxa descending from AncNode
#' @keywords internal
#' @export
tabulatetree = function(tree){

  # Tabulate branch information
  df = as.data.frame(tree$edge); colnames(df) = c("AncNode", "DescNode")
  df$Rbranch = 1:nrow(df) # add branch ID
  df$length = tree$edge.length # add branch length

  # Get descendants
  des = branchdefs(tree)

  # collapse to string
  des = unlist(lapply(des, function(x)paste0(sort(x), collapse = ",")))
  des = data.frame(DescNode = 1:length(des), des = des)

  # Add to table
  df = merge(df, des)

  # Re-order
  df = df[order(df$Rbranch),]

  # Now also add all branches descending to the table, too
  alldes = lapply(df$des, function(x)unlist(strsplit(x,split=",")))
  df$NoDes = unlist(lapply(alldes,length))
  df$desbranches = unlist(lapply(alldes,function(y)paste0(df$Rbranch[which(unlist(lapply(alldes,function(x)all(x %in% y))))], collapse = ",")))



  return(df)
}

#' @title Calculate path length across a tree.
#' @description Extract the root to tip distances for every terminal taxa (and optionally, internal node) in the tree.
#' @param tree A tree of class \code{phylo} used as the reference tree.
#' @param nodes Boolean operator (default = \code{FALSE}) that if \code{TRUE} will retain path length estimates for internal nodes as well as terminal taxa.
#' @return A named numeric vector of root-to-tip distances.
#' @importFrom ape dist.nodes Ntip Nnode
#' @export
PL = function(tree, nodes = F){
  pl = as.matrix(dist.nodes(tree))[,length(tree$tip)+1]
  names(pl) = c(rbind(tree$tip.label), (Ntip(tree)+1):(Ntip(tree)+Nnode(tree)))

  if(nodes == F)
    pl = pl[1:Ntip(tree)]

  return(pl)
}

#' @title Identify descendants of a node as named tips
#' @description Wrapper for phytools function getDescendants that retrieves tip labels rather than just node labels.
#' @param tree a phylogenetic tree as an object of class \code{phylo}.
#' @param node an integer specifying a node number in the tree.
#' @return The set of named tips descended from node in a vector.
#' @importFrom phytools getDescendants
#' @export
tipDescendants = function(tree, node){
  des = getDescendants(tree, node)
  des = des[des <= Ntip(tree)]
  des = tree$tip.label[des]
  return(des)
}


#' @title Compare branch lengths between two trees of identical topology
#' @description A function that extracts branch lengths between two trees of identical topology but where (for whatever reason) the edge matrices of the two trees do not match. This problem is encountered often when calculating stretched trees from the variable rates post-processor and comparing them to  means calculated from the raw trees output by the variable rates model.
#' @param tree1 a phylogenetic tree as an object of class \code{phylo}.
#' @param tree2 a phylogenetic tree as an object of class \code{phylo}.
#' @return The set of named tips descended from node in a vector.
#' @importFrom adephylo listTips
#' @importFrom graphics abline
#' @export
compareBLS = function(tree1, tree2){

  # Check
  if(!all(tree1$tip.label %in% tree2$tip.label)) stop("The two supplied trees do not contain the same taxa. ")

  # Extract descendants for each edge
  t1_descendants <- listTips(tree1)
  t2_descendants <- listTips(tree2)

  # Create unique identifiers for each branch
  t1_identifiers <- sapply(t1_descendants, function(x)paste0(sort(x),collapse = ","))
  t2_identifiers <- sapply(t2_descendants, function(x)paste0(sort(x),collapse = ","))

  # Data frame
  t1 = data.frame(id = t1_identifiers, t1node = (Ntip(tree1)+1):(Ntip(tree1)+Nnode(tree1)))
  t2 = data.frame(id = t2_identifiers, t2node = (Ntip(tree2)+1):(Ntip(tree2)+Nnode(tree2)))

  # Add branch lengths
  t1$t1brlen = tree1$edge.length[match(t1$t1node, tree1$edge[,2])]
  t2$t2brlen = tree2$edge.length[match(t2$t2node, tree2$edge[,2])]

  # merge together
  compared = merge(t1,t2, by = "id")

  # Plot the corresponding branch lengths
  plot(as.numeric(compared$t1brlen), as.numeric(compared$t2brlen), xlab = "Tree 1 Branch Lengths", ylab = "Tree 2 Branch Lengths", main = "Branch Length Comparison")
  abline(0, 1, col = "red")  # Optional: Add a y = x reference line

  return(compared)
}


#' @title Identify tippytomies
#' @description A function that identifies (and optionally deletes) 'tippytomies' from a tree. A tippytomy refers to the situation in which a pair (or group) of taxa have very short terminal branches resulting in what is essentially a bifurcation or polytomy at the tips of the tree. This can cause numerical issues in several programs. For instance, if two sister taxa have very short branches but very different trait values, the variable rates model will capitalize on this by maximizing the rate of evolution on these branches.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param cutoff The branch length below which taxa will be identified / removed.
#' @param action One of two character inputs: "identify" or "delete". See Returns for Details.
#' @return If "identify", the function will return a list of taxa that belong to each tippytomy. If "delete", the list is returned along with a tree which has deleted all but one (the first) member of each tippytomy. This tree can be used in subsequent comparative analyses.
#' @importFrom ape drop.tip
#' @export
tippytomies = function(tre, cutoff = 0.1, action = "identify"){
  # Identify branches that are less than the cut-off length
  shortbranches = which(tre$edge.length < cutoff)

  # Identify which of these branches are terminal branches
  shortterminals = shortbranches[which(tre$edge[shortbranches,2] <= Ntip(tre))]

  # Sort these numerically
  shortterminal_taxa = sort(tre$edge[shortterminals,2])

  # Identify groups of related branches
  res = split(shortterminal_taxa, cumsum(c(1,diff(shortterminal_taxa) !=1)))

  if(action == "identify") {
    res = lapply(res, function(x)tre$tip.label[x])
    return(res)
  }

  if(action == "delete"){
    # keep only the first in each cluster
    delete = unlist(lapply(res,function(x)tre$tip.label[x[2:length(x)]]))
    tre = drop.tip(tre, delete)
    tre$node.label = NULL
    res = lapply(res, function(x)tre$tip.label[x])
    return(list(tre = tre, tippytomies = res))
  }
}
