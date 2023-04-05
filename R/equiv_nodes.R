#' @title Get equivalences of nodes between phylogenies
#'
#' @description For a given set of trees of close construction (i.e., tips in common), this function returns the equivalences in node labels (or numbers if no node labels are provided) for each node given the subtaxa of each
#'
#' @usage equiv.nodes(trees,nodes,equitips=NULL,out=c("matrix","vector"),out.collapse="/")
#'
#' @param trees Phylogenetic trees to consider
#' @param nodes Optional. The nodes to be compared; without specifying it, all nodes of all trees are compared. Can be the node number (from 1 to the number of nodes, or from the number of tips +1 to the number of tips + number of nodes)
#' @param equitips Optional. If there are "equivalent tips" between phylogenies (i.e., same tip but not same label for that tip), they should be provided as a list of vectors, each containing the equivalent labels (the first one being to be used)
#' @param out The type of output. Can be a matrix (if \code{out="matrix"}, the default) with nodes as columns, phylogenies as rows, and node labels (or numbers) of the corresponding node for the given "global" node and phylogeny. Can be a vector (if \code{out="vector"}) collating node labels (or numbers).
#' @param out.collapse The way to collate node labels (or numbers) if a vector is desired as output. Default is set to \code{"/"}
#'
#' @import ape
#' @importFrom stats setNames
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get three random trees of 20, 19, and 21 tips (hence with some tips not present in all trees)
#' ## Start by having the larger tree, of 21 tips
#' set.seed(1)
#' tree<-rtree(21)
#' plot(tree)
#' axisPhylo()
#' ## Removing taxon #10 to get the first other tree
#' arrows(nodeHeights(tree)[tree$edge[,2]==10,2]+0.3,10,nodeHeights(tree)[tree$edge[,2]==10,2]+0.13,10,length=0.1,lwd=3,col="red")
#' ## Removing the taxa #5 and #14 for the second other tree
#' for(i in c(5,14)){
#'   arrows(nodeHeights(tree)[tree$edge[,2]==i,2]+0.3,i,nodeHeights(tree)[tree$edge[,2]==i,2]+0.13,i,length=0.1,lwd=3,col="green3")
#' }
#' trees<-c(drop.tip(tree,10),drop.tip(tree,c(5,14)),tree)
#'
#' # Getting node equivalences
#' equiv.nodes(trees)
#' # Checking node names and equivalence
#' par(mfrow=c(1,3),mar=rep(0,4))
#' for(i in 1:3){
#'   plot(trees[[i]],show.tip.label = FALSE)
#'   nodelabels()
#'   tiplabels()
#' }
#' @export

equiv.nodes<-function(trees,nodes,equitips=NULL,out=c("matrix","vector"),out.collapse="/"){
  if(missing(nodes)){
    nodes<-c(1:max(unlist(lapply(trees,Nnode))))
  }
  else if(is.character(nodes)){
    nodes<-c(1:length(which(nodes%in%unique(unlist(lapply(trees,function(x){trees$node.labels}))))))
  }

  if(length(out)>1){
    out="matrix"
  }

  all_combs<-lapply(trees,function(y){unlist(lapply(nodes,function(x){
    ST<-supporting.tips(y,x+Ntip(y),"comb")
    ST<-setNames(ST,rep(x,length(ST)))
  }))})

  if(!is.null(equitips)){
    for(i in 1:length(all_combs)){
      for(j in 1:length(all_combs[[i]])){
        for(k in 1:length(equitips)){
          tax<-strsplit(all_combs[[i]][j]," - ")[[1]]
          if(any(tax%in%equitips[[k]][-1])){
            tax[tax%in%equitips[[k]][-1]]<-equitips[[k]][1]
            all_combs[[i]][j]<-paste0(tax,collapse=" - ")
          }
        }
      }
    }
  }

  for(i in 1:length(trees)){
    if(is.null(trees[[i]]$node.label)){
      trees[[i]]$node.label<-as.numeric(c(1:Nnode(trees[[i]]))+Ntip(trees[[i]]))
    }
  }

  equiv<-matrix(ncol=length(nodes),nrow=length(trees))
  ref<-which.max(unlist(lapply(trees,Nnode)))
  if(Nnode(trees[[ref]])<length(nodes)){
    equiv[ref,]<-c(trees[[ref]]$node.label,rep(NA,length(nodes)-Nnode(trees[[ref]])))
  }
  else{
    equiv[ref,]<-trees[[ref]]$node.label
  }
  for(i in c(1:length(trees))[-ref]){
    for(j in 1:length(nodes)){
      local_equiv<-which(all_combs[[i]]%in%all_combs[[ref]][names(all_combs[[ref]])==as.character(j)])
      if(length(local_equiv)==0){
        equiv[i,j]<-NA
      }
      else{
        equiv[i,j]<-trees[[i]]$node.label[unique(as.numeric(names(all_combs[[i]][local_equiv])))]
      }
    }
  }

  if(out=="vector"){
    equiv<-apply(equiv,2,function(x){paste0(x[!is.na(x)],collapse=out.collapse)})
  }
  return(equiv)
}
