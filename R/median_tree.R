#' @title Get the tree with median taxa ages among a list of trees
#'
#' @description For a given list of trees of same structure (same tips, same branching pattern) but of various ages of taxa, this function computes (and returns) a single tree with median node ages.
#' It also names the nodes by either taking the already provided node labels or by creating new ones.
#'
#' @param trees The list of phylogenetic trees to consider
#' @param node.labelling Logical. Whether to set node labels or not.
#' @param node.label.opt A list of arguments to be passed to the \link[ape]{makeNodeLabel} function
#'
#' @import ape
#' @importFrom phytools nodeHeights
#' @importFrom stats median
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get a random tree of 20 tips
#' set.seed(1)
#' tree<-rtree(20)
#' tree<-ULT::add.ages.phylo(tree,tip.ages=rep(0,Ntip(tree)))
#' # Create 50 trees with some variation in its taxa ages (except the root), between the actual age and parent node age
#' nH<-max(nodeHeights(tree))-nodeHeights(tree) # Get ages (from the root) for each taxon
#' ages_bounds<-t(mapply(function(x){
#'   if(x==(Ntip(tree)+1)){
#'     rep(max(nH),2)
#'   }
#'   else{
#'     c(nH[tree$edge[,2]==x,2],nH[tree$edge[,2]==x,1])
#'   }
#' },c(1:(Ntip(tree)+Nnode(tree))))) # Get bounds of variation for each taxon
#' trees<-list()
#' for(i in 1:50){
#'   set.seed(i+1)
#'   trees[[i]]<-tree
#'   new_nH<-nH
#'   for(j in 1:(Ntip(tree)+Nnode(tree))){
#'     new_age<-runif(1,ages_bounds[j,1],ages_bounds[j,2])
#'     for(k in 1:2){
#'       new_nH[which(tree$edge[,k]==j),k]<-new_age
#'     }
#'   }
#'   new_nH<-max(new_nH)-new_nH
#'   trees[[i]]$edge.length<-apply(new_nH,1,diff)
#' }
#' # Plot four of them to just check the discrepancies in taxa ages
#' par(mfrow=c(2,2),mar=c(2.1,0,0,0))
#' for(i in 1:4){
#'   plot(trees[[i]])
#'   axisPhylo()
#' }
#' # Get tree with median taxa ages
#' medtree<-median.tree(trees)
#'
#' # Plot median tree
#' par(mfrow=c(1,1))
#' plot(medtree,show.tip.label = FALSE,x.lim=c(0,4.4))
#' for(i in 1:Ntip(medtree)){
#'   text(4.4,i,medtree$tip.label[i],adj=c(0,0.5),font=3)
#' }
#' axisPhylo()
#' # Add bars for the age variation of each taxon (light gray for tips, darker gray for nodes)
#' phylo.coords<-get("last_plot.phylo", envir = ape::.PlotPhyloEnv)[c("root.time","yy")]
#' for(i in 1:nrow(ages_bounds)){
#'   if(i==(Ntip(medtree)+1)){next}
#'   segments(phylo.coords$root.time-ages_bounds[i,1],phylo.coords$yy[i],phylo.coords$root.time-ages_bounds[i,2],phylo.coords$yy[i],col=rgb(t(rep(ifelse(i<Ntip(medtree),0.5,0),3)),alpha=0.25),lwd=10)
#' }
#' # Add original tree as shaded edges, to check whether bounds have been respected well
#' par(new=TRUE)
#' plot(tree,edge.col=rgb(1,0,0,0.25),show.tip.label=FALSE,x.lim=c(0,4.4))
#' par(mar=c(5,4,4,2)+0.1)
#'
#' @export median.tree

median.tree<-function(trees,node.labelling=TRUE,node.label.opt=NULL){
  med_tree<-trees[[1]]
  med_NH<-matrix(apply(mapply(nodeHeights,trees),1,median),ncol=2,nrow=nrow(med_tree$edge),byrow=FALSE)
  med_tree$edge.length<-apply(med_NH,1,diff)
  if(node.labelling&all(is.null(med_tree$node.labels))){
    mNL.args<-node.label.opt[names(node.label.opt)%in%c("method","prefix","nodeList")]
    mNL.args<-c(list("phy"=med_tree),mNL.args)
    med_tree<-do.call("makeNodeLabel",mNL.args)
  }
  return(med_tree)
}
