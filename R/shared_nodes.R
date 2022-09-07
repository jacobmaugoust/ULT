#' @title Getting shared nodes between several phylogenies
#'
#' @description This function returns, for each phylogeny, a list of logical vectors for the nodes being shared between different topologies
#'
#' @param trees The initial variable
#' @param own Optional. A list giving a subset of retained taxa for each phylogeny, if one does not want to consider all phylogenetic tips
#' @param shared Optional. A vector giving the common taxa between phylogenies to be retained. Automatically computed if not given as being all shared tips.
#'
#' @import ape
#'
#' @examples
#' set.seed(1)
#' orig.tree<-rtree(50)
#' if(!require(foreach)){
#'   subsets<-list()
#'   trees<-list()
#'   for(i in 1:4){
#'     subsets[[i]]<-sample(orig.tree$tip.label,30,FALSE)
#'     trees[[i]]<-keep.tip(orig.tree,subsets[[i]])
#'   }
#' } else {
#'   subsets<-foreach(i=1:4)%do%{sample(orig.tree$tip.label,30,FALSE)}
#'   trees<-foreach(i=1:4)%do%{keep.tip(orig.tree,subsets[[i]])}}
#' Reduce(intersect,subsets) # Common species - only 7 species, that's convenient to see common nodes
#' par(mfrow=c(2,2),mar=c(3.1,2.1,2.1,0.1))
#' for(i in 1:4){ # Let's now plot the trees, and highlight the common species, the edges leading to them, and the nodes involved
#'   if(require(RRphylo)){
#'     inv.nodes<-unlist(sapply(Reduce(intersect,subsets),function(x){RRphylo::getMommy(trees[[i]],which(trees[[i]]$tip.label%in%x))}))
#'     plot(trees[[i]],show.tip.label=FALSE,edge.color=binarize(trees[[i]]$edge[,2],zero=c(unique(inv.nodes),which(trees[[i]]$tip.label%in%Reduce(intersect,subsets))),output=c("red","black")))
#'   } else {plot(trees[[i]],show.tip.label=FALSE)}
#'   tiplabels(trees[[i]]$tip.label[trees[[i]]$tip.label%in%Reduce(intersect,subsets)],
#'             which(trees[[i]]$tip.label%in%Reduce(intersect,subsets)),adj=-0.1,col="red",frame="none",font=3,cex=0.7)
#'   tiplabels(trees[[i]]$tip.label[!trees[[i]]$tip.label%in%Reduce(intersect,subsets)],
#'             which(!trees[[i]]$tip.label%in%Reduce(intersect,subsets)),adj=-0.1,col="black",frame="none",font=3,cex=0.7)
#'   if(require(RRphylo)){
#'     inv.dic.nodes<-unique(inv.nodes)[sapply(unique(inv.nodes),function(x){length(which(trees[[i]]$edge[trees[[i]]$edge[,1]==x,2]%in%c(unique(inv.nodes),which(trees[[i]]$tip.label%in%Reduce(intersect,subsets)))))})>1]
#'     nodelabels(as.character((Ntip(trees[[i]])+1):(Ntip(trees[[i]])+Nnode(trees[[i]])))[inv.dic.nodes-Ntip(trees[[i]])],
#'              c((Ntip(trees[[i]])+1):(Ntip(trees[[i]])+Nnode(trees[[i]])))[inv.dic.nodes-Ntip(trees[[i]])],col="red",cex=0.7,adj=-0.5,frame="none")
#'     nodelabels(as.character((Ntip(trees[[i]])+1):(Ntip(trees[[i]])+Nnode(trees[[i]])))[-(inv.dic.nodes-Ntip(trees[[i]]))],
#'              c((Ntip(trees[[i]])+1):(Ntip(trees[[i]])+Nnode(trees[[i]])))[-(inv.dic.nodes-Ntip(trees[[i]]))],col="black",cex=0.7,adj=-0.5,frame="none")
#'   } else{nodelabels(frame="none",col="black",cex=0.7,adj=-0.5)}
#' }
#' shared.nodes(trees) # Here is the direct output of the function, with logicals for nodes of each tree
#' for(i in 1:4){print(as.character(c((Ntip(trees[[i]])+1):(Ntip(trees[[i]])+Nnode(trees[[i]]))))[shared.nodes(trees)[[i]]])} # Here are the node numbers retained for each tree (compare with the plot)
#'
#' @export

shared.nodes<-function(trees,own=NULL,shared=NULL){
  if(!is.null(own)){
    for(i in 1:length(trees)){
      trees[[i]]<-keep.tip(trees[[i]],own[[i]])
    }
  }
  for(i in 1:length(trees)){
    trees[[i]]$node.label<-as.character((Ntip(trees[[i]])+1):(Nnode(trees[[i]])+Ntip(trees[[i]])))
  }
  if(is.null(shared)){
    shared<-do.call("Reduce",list("intersect",lapply(trees,function(x){x$tip.label})))
  }
  final<-list()
  for(i in 1:length(trees)){
    final[[i]]<-trees[[i]]$node.label%in%keep.tip(trees[[i]],shared)$node.label
  }
  return(final)
}
