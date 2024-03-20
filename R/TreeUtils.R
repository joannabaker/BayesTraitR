#' @title Calculate path length across a tree.
#' @description Extract the root to tip distances for every terminal taxa (and optionally, internal node) in the tree.
#' @param phy A tree of class("Phylo") used as the reference tree.
#' @param nodes Boolean operator (default value F) that if TRUE will retain path length estimates for internal nodes as well as terminal taxa.
#' @return A named numeric vector of root-to-tip distances.
#' @importFrom ape dist.nodes Ntip Nnode
#' @export
PL = function(phy, nodes = F){
  pl = as.matrix(dist.nodes(phy))[,length(phy$tip)+1]
  names(pl) = c(rbind(phy$tip.label), (Ntip(phy)+1):(Ntip(phy)+Nnode(phy)))

  if(nodes == F)
    pl = pl[1:Ntip(phy)]

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

