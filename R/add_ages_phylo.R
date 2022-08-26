#' @title Adding taxa and/or node ages to a phylogeny
#'
#' @description This function transforms a \code{phylo} object to add to it taxa ages at the nodes and/or at the terminal tips
#'
#' @param tree A phylogenetic tree.
#' @param tip.ages Optional. A numeric vector with ages of the terminal tips. Must be of same length than the number of tree tips. By default, values are considered to be in the order of the tree tips; this is bypassed if the vector is named with the tip labels.
#' @param node.ages Optional. A numeric vector with ages of the internal nodes. By default, values are considered to be in the order of the tree nodes; this is bypassed if the tree has node labels and if this vector is named with these node labels.
#' @param plot Optional. Turned to \code{TRUE} by default, meaning that the phylogeny is plotted at the end of the execution of the function. Turn to \code{FALSE} if not desired.
#'
#' @importFrom phytools nodeHeights
#' @import ape
#'
#' @examples
#' require(ape)
#' # Create an ultrametric tree
#' tree<-rtree(20)
#' tree<-tree[-which(names(tree)=="edge.length")]
#' class(tree)<-"phylo"
#' # Create random tip ages
#' tip.ages<-runif(Ntip(tree),0,5)
#' # Create random node ages, sorting them so that the more basal a node is, the oldest it is (and therefore avoid negative branch lengths)
#' node.ages<-sort(runif(Nnode(tree),6,10),decreasing = TRUE)
#' # Add tip and node ages to the tree
#' newtree<-add.ages.phylo(tree,tip.ages,node.ages)
#'
#' @export

add.ages.phylo<-function(tree,tip.ages=NULL,node.ages=NULL,plot=TRUE){
  if(!is.null(tip.ages)){
    if(all(is.null(names(tip.ages)))){
      names(tip.ages)<-tree$tip.label
    }
  }

  is.node.label<-TRUE
  if(!"node.label"%in%names(tree)){
    is.node.label<-FALSE
    tree<-makeNodeLabel(tree,prefix="N")
  }

  if(!is.null(node.ages)){
    if(all(is.null(names(node.ages)))){
      names(node.ages)<-tree$node.label
    }
  }

  named.edge<-tree$edge
  for(i in 1:nrow(tree$edge)){
    for(j in 1:ncol(tree$edge)){
      named.edge[i,j]<-c(tree$tip.label,tree$node.label)[tree$edge[i,j]]
    }
  }

  nodes<-list()
  for(i in 1:Nnode(tree)){
    rev_i<-Nnode(tree)+1-i
    nodes[[rev_i]]<-named.edge[which(named.edge[,1]==tree$node.label[i]),2]
    term.br<-which(!is.na(suppressWarnings(as.numeric(nodes[[rev_i]]))))
    if(length(term.br)>0){
      for(j in 1:length(term.br)){
        nodes[[rev_i]][term.br[j]]<-tree$tip.label[as.numeric(nodes[[rev_i]][term.br[j]])]
      }
    }
  }
  names(nodes)<-tree$node.label
  anc.nodes<-nodes
  for(i in 1:Nnode(tree)){
    int.br<-which(anc.nodes[[i]]%in%tree$node.label)
    if(length(int.br)>0){
      for(j in 1:length(int.br)){
        nodes[[i]][int.br[j]]<-names(anc.nodes)[which(mapply(function(x){all(x==c(tree$tip.label,tree$node.label)[tree$edge[tree$edge[,1]==Ntip(tree)+which(tree$node.label==anc.nodes[[i]][int.br[j]]),2]])},anc.nodes))]
      }
    }
  }
  names(nodes)<-tree$node.label
  nodes<-lapply(nodes,function(x){rev(x)})

  NT_arguments<-list(taxa=tree$tip.label,
                     nodes=nodes,
                     plot=plot)
  if(!is.null(tip.ages)){
    NT_arguments<-c(NT_arguments,list(age_taxa=tip.ages[match(tree$tip.label,names(tip.ages))]))
  }
  if(!is.null(node.ages)){
    NT_arguments<-c(NT_arguments,list(age_nodes=node.ages[match(tree$node.label,names(node.ages))]))
  }
  if(is.node.label){
    NT_arguments<-c(NT_arguments,node.labels=tree$node.label)
  }
  else{
    NT_arguments<-c(NT_arguments,node.labels=FALSE)
  }
  new.tree<-do.call("create.tree",NT_arguments)
  if(plot){
    axisPhylo()
  }
  return(new.tree)
}

