#' @title Substitute a taxon of a given phylogenetic tree by another tree
#'
#' @description The goal of this function is to enhance manipualtion between trees, especially for constructing composite trees, by allowing one to substitute a taxon (tip or subtree) of a tree by another tree.
#'
#' @details
#' This function replaces a taxon of a given "background" tree, be it a single tip or a node (i.e., part of the tree), by another "foreground" tree.
#' The age of the node of the "foreground" tree is not necessarily the same than the age of the taxon to replace in the "background" tree; there are therefore options allowing for specifying its age.
#'
#' @param back.tree The tree on which a taxon is substituted by another tree, named "background" tree because it is the general support/scaffold/backbone.
#' @param fore.tree The tree to put in the "background" tree, named "foreground" tree because it is the one to put on the other tree.
#' @param where The taxon of \code{back.tree} one wants to replace by the \code{fore.tree}. It can be the number or the name of the taxon (tip or node), but it can also be a vector of some tips (allowing for specifying a node without having to search for it and without having to quote every of its descendants).
#' @param drop.where Optional. The user can discard the given taxon where to put the \code{fore.tree} (default condition, set to \code{TRUE}), or it can keep it if the user wants to place the \code{fore.tree} on the edge leading to that taxon (if set to \code{FALSE}).
#' @param poly.where Optional. The user can chose whether the \code{fore.tree} has to be on a distinct edge (default condition, set to \code{FALSE}) or if it has to be in a polytomy with the specified taxon if the taxon is a node (if set to \code{TRUE}).
#' @param stem.edge.length Optional. The length of the edge leading to the \code{fore.tree} in the new tree. Complementary with the \code{node.age} parameter; please do not specify both.
#' @param node.age Optional. The age of the 'root' node of the \code{fore.tree} is the new tree. Must be lower than the age of the previous node (if not, a warning message with the age of that previous node is displayed). Complementary with the \code{stem.edge.length} parameter; please do not specify both.
#'
#' @return A tree of class \code{phylo} with both the \code{back.tree} and the \code{fore.tree} merged together.
#'
#' @examples
#' back.tree<-read.tree(text="((((A:1,B:2):0.5,(C:3,(D:4,E:1):2):5):1,(F:3,(G:2.5,(H:2.5,(I:1.5,J:1):1):1):2):1):4,K:10);")
#' plot(ref_tree)
#' nodelabels()
#'
#' fore.tree<-read.tree(text="((x:0.5,y:1.5):5,z:11);")
#' plot(toswap_tree)
#'
#' where<-c("C","E")
#'
#' plot(subs.tree(back.tree,fore.tree,where,node.age=8))
#' plot(subs.tree(back.tree,fore.tree,where,stem.edge.length=0.1))
#' plot(subs.tree(back.tree,fore.tree,where))
#'
#' @import ape
#' @importFrom gtools ask
#' @importFrom phytools getDescendants nodeHeights

subs.tree<-function(back.tree,fore.tree,where,drop.where=TRUE,poly.where=FALSE,stem.edge.length=NULL,node.age=NULL){
  if(missing(back.tree)|missing(fore.tree)|missing(where)){
    error("Please specifiy 'back' and 'fore' trees and a point on the 'back' tree to insert the 'fore' one")
  }

  where<-ifelse(length(where)>1,ULT::getClade(back.tree,where),where)
  where<-ifelse(is.character(where),ifelse(where%in%back.tree$tip.label,which(back.tree$tip.label==where),which(back.tree$where.label==where)),where)
  if(where>(Ntip(back.tree)+Nnode(back.tree))){
    error("Number provided is neither that of a tip nor that of a node")
  }

  if(!missing(stem.edge.length)&!missing(node.age)){
    choice<-ask("Both a length for the edge of the stem and an age for the node have been provided, which one do you want to keep? (stem/node)")
    if(choice=="stem"){
      node.age<-NULL
    }
    else if(choice=="node"){
      stem.edge.length<-NULL
    }
    else{
      rm(node.age, stem.edge.length)
    }
  }

  if(poly.where==TRUE&where<=Ntip(back.tree)){
    poly.where<-FALSE
  }

  if(where<=Ntip(back.tree)){
    tree<-bind.tree(x=back.tree,y=fore.tree,where=where)
  }
  else{
    tax_to_rm<-getDescendants(back.tree,where)
    edge_lines_to_rm<-which(back.tree$edge[,2]%in%tax_to_rm)
    edge_to_rm<-back.tree$edge[edge_lines_to_rm,]

    n_tips_diff<-Ntip(fore.tree)-length(tax_to_rm[tax_to_rm<=Ntip(back.tree)])
    n_nodes_diff<-Nnode(fore.tree)-length(tax_to_rm[tax_to_rm>Ntip(back.tree)])-1

    if(all(edge_lines_to_rm==1)){
      keep<-c(1,Ntip(back.tree)+1,Ntip(back.tree)+2,Ntip(back.tree)+Nnode(back.tree))
    }
    else{
      keep<-c(max(edge_to_rm[edge_to_rm<=Ntip(back.tree)])+1,
              min(edge_to_rm[edge_to_rm>Ntip(back.tree)]),
              max(edge_to_rm[edge_to_rm>Ntip(back.tree)])+1,
              Ntip(back.tree)+Nnode(back.tree))
    }

    renum_edge<-back.tree$edge
    renum_edge[back.tree$edge%in%c(keep[1]:keep[2])]<-back.tree$edge[back.tree$edge%in%c(keep[1]:keep[2])]+n_tips_diff
    renum_edge[back.tree$edge%in%c(keep[3]:keep[4])]<-back.tree$edge[back.tree$edge%in%c(keep[3]:keep[4])]+n_tips_diff+n_nodes_diff

    before<-renum_edge[min(1,min(edge_lines_to_rm)-1):min(edge_lines_to_rm)-1,]

    between<-fore.tree$edge
    between[between>Ntip(fore.tree)]<-between[between>Ntip(fore.tree)]-Ntip(fore.tree)-1+max(before)
    between[between<=Ntip(fore.tree)]<-between[between<=Ntip(fore.tree)]+max(0,before[before<c(before)[1]])

    after<-renum_edge[min(dim(renum_edge)[1],max(edge_lines_to_rm)+1):max(dim(renum_edge)[1],edge_lines_to_rm),]

    edge<-rbind(before,between,after)

    tip.label<-c(back.tree$tip.label[min(1,(min(tax_to_rm)-1)):(min(tax_to_rm)-1)],
                 fore.tree$tip.label,
                 back.tree$tip.label[keep[1]:Ntip(back.tree)])

    Nnode<-max(edge)-(edge[1,1]-1)

    tree<-list("edge"=edge,"Nnode"=Nnode,"tip.label"=tip.label)

    if(all(any(names(back.tree)=="edge.length"),any(names(fore.tree)=="edge.length"))){
      edge.length<-c(back.tree$edge.length[min(1,min(edge_lines_to_rm)-1):min(edge_lines_to_rm)-1],
                     fore.tree$edge.length,
                     back.tree$edge.length[min(dim(renum_edge)[1],max(edge_lines_to_rm)+1):max(dim(renum_edge)[1],edge_lines_to_rm)])
      tree<-c(tree,list("edge.length"=edge.length))
    }
    attr(tree,"class")<-"phylo"
    attr(tree,"order")<-"cladewise"
  }

  if(any(names(tree)=="edge.length")){
    stem<-which(back.tree$edge[,2]==where)
    new.edge.length<-tree$edge.length

    if(!is.null(stem.edge.length)){
      new.edge.length[stem]<-stem.edge.length
    }
    if(!is.null(node.age)){
      prev.stem.age<-(max(nodeHeights(back.tree))-nodeHeights(back.tree)[which(back.tree$edge[,2]==back.tree$edge[stem,1]),2])
      if(node.age>prev.stem.age){
        warning(paste0("Node age provided exceeds that of the parent node (i.e., age value: ",prev.stem.age,"); node age reset"))
        node.age<-(max(nodeHeights(back.tree))-nodeHeights(back.tree)[stem,2])
      }
      new.edge.length[stem]<-new.edge.length[stem]+(max(nodeHeights(back.tree))-nodeHeights(back.tree)[stem,2])-node.age
    }

    tree$edge.length<-new.edge.length
  }
  return(tree)
}

temp_fun<-function(x,y,z,aaa){
  temp<-match.call()
  return(as.list(temp))
}
test<-temp_fun(1,2,3,4)
test
?match.call
