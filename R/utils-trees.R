#' @title Get all descendants of all nodes in a tree
#' @description Retrieves the descendants of every node (including terminals) in a tree.
#' @param tree a phylogenetic tree as an object of class "phylo".
#' @return A list of vectors defining the taxa descending from every node of a tree.
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
#' @param tree a phylogenetic tree as an object of class "phylo".
#' @return A data.frame with n rows equal to the number of branches in the tree
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
  des = unlist(lapply(des, function(x)paste0(sort(x), collapse = ",")))
  des = data.frame(DescNode = 1:length(des), des = des)

  # Add to table
  df = merge(df, des)

  # Re-order
  df = df[order(df$Rbranch),]

  # Now also add all branches descending to the table, too
  alldes = lapply(df$des, function(x)unlist(strsplit(x,split=",")))
  df$desbranches = unlist(lapply(alldes,function(y)paste0(df$Rbranch[which(unlist(lapply(alldes,function(x)all(x %in% y))))], collapse = ",")))

  return(df)
}

#' @title Calculate path length across a tree.
#' @description Extract the root to tip distances for every terminal taxa (and optionally, internal node) in the tree.
#' @param tree A tree of class("Phylo") used as the reference tree.
#' @param nodes Boolean operator (default value F) that if TRUE will retain path length estimates for internal nodes as well as terminal taxa.
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
#' @param tree a phylogenetic tree as an object of class "phylo".
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
