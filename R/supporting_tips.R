#' @title Get all pairs of tips defining a node
#'
#' @description For a given tree and a given node of this tree, this function outputs all pairs (or more if polytomies) of tips whose MRCA (most recent common ancestor) is that node of that tree
#'
#' @param tree Phylogenetic tree to consider
#' @param node Node to consider. Can be the node number (from 1 to the number of nodes, or from the number of tips +1 to the number of tips + number of nodes)
#' @param out The type of output. Can be a list (if \code{out="list"}, the default) containing vectors of tip labels; can be a vector (if \code{out="comb"}) collating tip labels
#' @param out.collapse The way to collate tip labels if a vector is desired as output. Default is set to \code{" - "}
#' @param mandatory Logical, whether "mandatory tips" (tips always found to structure the node, likely one of the direct daugher taxa of the given node) should be returned. Default to \code{FALSE}
#'
#' @import ape
#'
#' @examples
#' require(ape)
#' # Get a random tree of 20 tips
#' set.seed(1)
#' tree<-rtree(20)
#' plot(tree)
#' nodelabels()
#' # To see which node is "structuring" the node 27
#' supporting.tips(tree,27)
#' # Same with the vectorized output
#' supporting.tips(tree,27,"comb")
#' # One of the direct daughter taxa of the node 27 is the taxon t17, so it is a "mandatory" tip (i.e., always in the combinations)
#' supporting.tips(tree,27,mandatory=TRUE)
#' @export

supporting.tips<-function(tree,node,out=c("list","comb"),out.collapse=" - ",mandatory=FALSE){
  if(all(is.null(tree$tip.label))){
    stop("Please provide a tree with tip labels")
  }
  if(!is.numeric(node)){
    if(all(is.null(tree$node.label))){
      stop("Please provide either a node number or a tree with node labels corresponding to the inputted node name")
    }
    else{
      node<-Ntip(tree)+which(tree$node.label==node)
    }
  }
  else if(node<Ntip(tree)){
    node<-node+Ntip(tree)
  }
  if(length(out)==2){
    out<-"list"
  }

  tips_combs<-expand.grid(tree$tip.label,tree$tip.label)
  tips_combs<-unique(unlist(lapply(apply(tips_combs,1,function(x){if(x[1]==x[2]){NA}else{sort(x)}}),function(x){if(!any(is.na(x))){paste0(x,collapse=out.collapse)}})))
  good_combs<-tips_combs[sapply(tips_combs,function(x){node==getMRCA(tree,strsplit(x,out.collapse)[[1]])})]
  if(out=="list"){
    good_combs<-lapply(good_combs,function(x){strsplit(x,out.collapse)[[1]]})
  }
  ret<-good_combs
  if(mandatory){
    temp<-table(unlist(good_combs))
    if(any(temp==length(good_combs))){
      mand<-names(temp[temp==length(good_combs)])
    }
    else{
      mand<-NA
    }
    ret<-c("Tip combinations"=list(good_combs),"Mandatory tips"=list(mand))
  }
  return(ret)
}
