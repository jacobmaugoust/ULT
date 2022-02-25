#' @title
#' Get the node of a given monophyletic group.
#'
#' @description
#' This function aims to output the monophyletic group (node number or name) that contains a given set of terminal taxa, and only them.
#'
#' @details
#' This function heavy relies on the function \code{getDescendants} from the package phytools, which give all phylogenetical units (i.e., tips and nodes) that are contained in a given clade.
#' By specifying only some of the descendants of a node, the function returns the apicalmost node uniting them (i.e., the crown group node).
#'
#' @param tree The phylogenetic tree to provide
#' @param taxa The taxa (tips only) forming the monophyletic group (i.e., a monophyletic group is at least composed of two tips).
#' @param output Optional. The format of the output (i.e., the node). Can be either the node number according to the \code{nodelabels} function (\code{output="number"}, the default), or its number across nodes (\code{output="node number"}), or its name if this has been specified in the tree (\code{output="name"}).
#'
#' @return
#' The number or name of the node containing the given tips.
#'
#' @importFrom phytools getDescendants
#' @import ape
#'
#' @examples
#' newick_tree<-c("((((A,B),C),(D,E)),F);") # Create a simple tree in NEWICk format
#' tree<-ape::read.tree(text = newick_tree) # Read it
#' plot(tree)
#' ape::nodelabels() # Plot it together with node numbers
#'
#' taxa<-c("A","B","C") # A monophyletic group
#' getClade(tree,taxa) # The actual 9 on the phylogeny
#' getClade(tree,taxa,output="node number") # The node number three going from the root
#'
#' tree2<-ape::makeNodeLabel(tree,method="number",prefix="N") # Adding names to the nodes
#' plot(tree2)
#' ape::nodelabels(tree2$node.label)
#' getClade(tree2,taxa,output="name")
#'
#' @export

getClade<-function(tree,taxa,output="number"){
  if(missing(tree)){
    stop("Please provide a phylogenetic tree in the tree format")
  }
  if(missing(taxa)){
    stop("No taxa provided")
  }
  if(is.numeric(taxa)){
    if(any(!taxa%in%c(1:Ntip(tree)))){
      stop("You provided taxa number that are exceeding the number of tips; please provide only tip numbers or names")
    }
    taxa<-tree$tip.label[taxa]
  }
  else if (is.character(taxa)){
    if(all(taxa%in%tree$tip.label)==FALSE){
      taxa<-taxa[taxa%in%tree$tip.label]
    }
  }
  if(output=="name"&exists('node.label',where=tree)==FALSE){
    stop("No node names provided in the tree but a named node has been required; please change the output option or give names to the tree nodes")
  }

  n_taxa<-Ntip(tree)
  that_node<-c()
  for (i in 1:tree$Nnode){
    n_node<-n_taxa+i
    n_desc<-getDescendants(tree,n_node)
    desc<-tree$tip.label[n_desc[n_desc<=n_taxa]]
    if(all(taxa%in%desc)){
      that_node<-c(that_node,n_node)
    }
  }
  if(length(that_node)>1){
    l_desc<-c()
    for (i in 1:length(that_node)){
      l_desc<-c(l_desc,length(!getDescendants(tree,that_node[i])%in%taxa))
    }
    that_node<-that_node[which.min(l_desc)]
  }

  if(output=="number"){
    return(that_node)
  }
  if(output=="node number"){
    return(that_node-n_taxa)
  }
  if(output=="name"){
    return(tree$node.label[that_node-n_taxa])
  }
}
