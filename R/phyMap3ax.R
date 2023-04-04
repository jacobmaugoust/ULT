#' @title
#' Summarize three traits or axes as RGB colors and map them on a phylogeny
#'
#' @description
#' This functions first rescales three given traits are RGB color codes, then map them on a phylogeny, estimating or using given ancestral estimates.
#'
#' @details
#' The aim of this function is to provide a visually intuitive summary of the information given by three traits or three axes (for instance, the three first axes of a PCA).
#' To do this, the three variables are rescaled as 0-255 RGB values to summarize each point as a color code.
#' Then, the aim is to give the color information of each taxon on a phylogeny, to be able to visually compare the phylogenetic (i.e., according to branches) and the statistical (i.e., according to colours) clusterings.
#' To have a smooth effect, it is better to also have ancestral estimates. These can be given as parameters of the function, or they can be estimated using different methods.
#' \cr
#' The method for reconstructing ancestral states (if not provided) is by default the phylogenetic ridge regression (see \code{\link[RRphylo]{RRphylo}}).
#' Other methods are available, especially those relying on maximum lilkelihood estimations, such as \code{\link[ape]{ace}} (ML estimation under Brownian motion model),
#' \code{\link[phytools]{fastAnc}} (ML estimation using Felsenstein's contrasts),
#' \code{\link[phytools]{anc.ML}} (ML estimation under various models such as Brownian motion, Ornstein-Uhlenbeck process, "early burst" model...),
#' \code{\link[phytools]{anc.trend}} (ML estimation under Brownian motion model with a trend),
#' \code{\link[phytools]{anc.Bayes}} (using Bayesian MCMC under adjustable Brownian motion model).
#'
#' @importFrom stats setNames
#'
#' @examples
#' # Simulate three traits, then summarize and map then on a phylogeny
#'
#' n<-50 # Number of tips
#' require(ape)
#' require(phytools)
#' tree<-rtree(n) # Generating random tree
#' a<-runif(n,min=10,max=50) # Generating three different trait within a given range
#' b<-runif(n,min=1,max=5)
#' c<-runif(n,min=10,max=20)
#' a<-sort(a) # Sorting values is actually adding a phylogenetic effect that will be visually from bottom to top of the phylogram; not sorting them allows for non-phylogenetic distribution
#' b<-sort(b)
#' # c<-sort(c)
#' names(a)<-names(b)<-names(c)<-tree$tip.label # Assigning tip names to data
#' data<-data.frame(a,b,c) # Generating dataset
#' method<-"fastAnc" # Choosing a method for reconstructing ancestral states
#'
#' phyMap3ax(tree,data,method) # And here is the phylogram with the branches coloured according to the values of the three traits
#'
#' # Simulate a general phylogenetic effect with some "outliers"
#'
#' require(RRphylo)
#' n<-50 # Number of tips
#' tree<-rtree(n) # Generating random tree
#' a<-setBM(tree,s2=2,a=10,type="brown") # Generating three different traits with phylogenetic effect under Brownian motion regime
#' b<-setBM(tree,s2=20,a=1,type="brown")
#' c<-setBM(tree,s2=0.2,a=0,type="brown")
#' out<-sample(1:n,round(n/10,0),replace=FALSE) # Randomly choose 10% of taxa that have weird values
#' a[out]<-runif(length(out),min=min(a)-(max(a)-min(a))*2,max=max(a)+(max(a)-min(a))*2) # Assign to the 'weird taxa' values within a larger range than the original trait range
#' b[out]<-runif(length(out),min=min(b)-(max(b)-min(b))*2,max=max(b)+(max(b)-min(b))*2)
#' c[out]<-runif(length(out),min=min(c)-(max(c)-min(c))*2,max=max(c)+(max(c)-min(c))*2)
#' names(a)<-names(b)<-names(c)<-tree$tip.label # Assigning tip names to data
#' data<-data.frame(a,b,c) # Generating dataset
#' method<-"fastAnc" # Choosing a method for reconstructing ancestral states
#'
#' phyMap3ax(tree,data,method) # Here is the plot of the three traits
#' if(require(ULT)){tiplabels(ULT::binarize(c(1:n)%in%out,zero=FALSE,one=TRUE,output=c("","< HERE")),adj=-0.5,frame="none",col="red",font=2)} # And here are the 'weird taxa'
#'
#' # Simulate a general phylogenetic effect with some convergent taxa
#'
#' require(RRphylo)
#' n<-50 # Number of tips
#' tree<-rtree(n) # Generating random tree
#' a<-setBM(tree,s2=2,a=10,type="brown") # Generating three different traits with phylogenetic effect under Brownian motion regime
#' b<-setBM(tree,s2=20,a=1,type="brown")
#' c<-setBM(tree,s2=0.2,a=0,type="brown")
#' conv<-sample(1:n,round(n/10,0),replace=FALSE) # Randomly choose 10% of taxa that converge
#' out_a<-sample(c(runif(1,min=min(a)-(max(a)-min(a))*2,max=min(a)),runif(1,min=max(a),max=max(a)+(max(a)-min(a))*2)),1) # Choosing an extreme (i.e., largely outside from the variation range) of the trait 'a'
#' a[conv]<-out_a+runif(length(conv),min=out_a-(max(a)-min(a))/2,max=out_a+(max(a)-min(a))/2) # Randomly assign to the convergent taxa values floating around the extreme value
#' out_b<-sample(c(runif(1,min=min(b)-(max(b)-min(b))*2,max=min(b)),runif(1,min=max(b),max=max(b)+(max(b)-min(b))*2)),1)
#' b[conv]<-out_b+runif(length(conv),min=out_b-(max(b)-min(b))/2,max=out_b+(max(b)-min(b))/2)
#' out_c<-sample(c(runif(1,min=min(c)-(max(c)-min(c))*2,max=min(c)),runif(1,min=max(c),max=max(c)+(max(c)-min(c))*2)),1)
#' c[conv]<-out_c+runif(length(conv),min=out_c-(max(c)-min(c))/2,max=out_c+(max(c)-min(c))/2)
#' names(a)<-names(b)<-names(c)<-tree$tip.label # Assigning tip names to data
#' data<-data.frame(a,b,c) # Generating dataset
#' method<-"fastAnc" # Choosing a method for reconstructing ancestral states
#'
#' phyMap3ax(tree,data,method) # Here is the plot of the three traits
#' if(require(ULT)){tiplabels(ULT::binarize(c(1:n)%in%conv,zero=FALSE,one=TRUE,output=c("","< HERE")),adj=-0.5,frame="none",col="red",font=2)} # And here are the convergent taxa
#'
#' @param tree The tree of format \code{phylo} to map the data on
#' @param data The data to map on the tree. Must be a data frame or a matrix with three columns.
#' @param method Optional. The method to use to reconstruct the ancestral states, if those are not provided. Available methods are: \code{RRphylo} (the default), \code{ace}, \code{anc.ML}, \code{anc.trend}, \code{anc.Bayes}, and \code{fastAnc}. See details for these methods.
#' @param method.opt Optional. The options to the method used if some customization is desired. Must be a list of named elements, and these names should be the options of the previously called method.
#' @param anc.states Optional. The ancestral states to the data. As for \code{data}, must be a three-columned data frame or matrix.
#' @param res Optional. The resolution of the coloring of the tree edges (i.e., the highest the resolution, the smoothest the color gradient). Default is 1000.
#' @param plot.opt Optional. Options to be passed for the plot. See \code{\link[phytools]{contMap}} and \code{\link[phytools]{plot.contMap}} help pages.
#'
#' @importFrom RRphylo RRphylo
#' @importFrom utils tail
#'
#' @export
phyMap3ax<-function(tree,data,method,method.opt,anc.states,res,plot.opt){
  col_from_three<-function(x){
    first<-matrix(ncol=3,nrow=dim(x)[1])
    for (i in 1:3){
      first[,i]<-scales::rescale(x[,i],to=c(0,1))
    }
    second<-apply(first,1,function(x){rgb(t(x))})
    return(second)
  }
  if(dim(data)[2]>3){
    data<-data[,1:3]
    warning("Provided dataset with more than three columns; only three first columns used")
  }

  if(missing(method)){
    if(missing(anc.states)){
      method<-"RRphylo"
    }
    else{
      method<-FALSE
    }
  }

  if(missing(method.opt)){method.opt<-NULL}

  if(missing(anc.states)){
    old_data<-data
    data<-matrix(nrow=Ntip(tree),ncol=3,NA)
    anc.states<-matrix(ncol=3,nrow=Nnode(tree),NA)
    for (i in 1:3){
      if(length(names(old_data[,i]))==0){
        if(all(rownames(old_data)%in%tree$tip.label)&length(rownames(old_data))==Ntip(tree)){
          x<-setNames(old_data[,i],rownames(old_data))
        }
        else{
          x<-setNames(old_data[,i],tree$tip.label)
        }
      }
      else{
        x<-old_data[,i]
      }
      x<-x[order(match(names(x),tree$tip.label))]
      data[,i]<-x
      basics<-list("tree"=tree,"x"=x)
      if(method=="RRphylo"){names(basics)[2]<-"y"}
      if(method=="ace"){names(basics)[1]<-"phy"}

      temp<-invisible(do.call(method,c(basics,method.opt)))

      if(method=="RRphylo"){
        anc.states[,i]<-temp$aces[,1]
      }
      if(method%in%c("anc.ML","anc.trend","ace")){
        anc.states[,i]<-temp$ace
      }
      if(method=="anc.Bayes"){
        anc.states[,i]<-unlist(tail(temp$mcmc[,(1:Nnode(tree)+2)],n=1))
      }
      if(method=="fastAnc"){
        anc.states[,i]<-temp
      }
    }
  }
  colnames(anc.states)<-colnames(data)
  rownames(anc.states)<-c((Ntip(tree)+1):(Ntip(tree)+Nnode(tree)))

  full_data<-rbind(data,anc.states)
  cols<-col_from_three(full_data)
  names(cols)<-row.names(full_data)
  ord_cols<-c(cols[as.numeric(t(tree$edge))])

  x<-setNames(1:Ntip(tree),tree$tip.label)
  y<-(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
  if(missing(res)){res<-1000}
  temp_CM<-contMap(tree,x=x,res=res,method="user",anc.states=y,plot=FALSE)

  nsteps_edges<-mapply(function(x){sum(table(names(x)))},temp_CM$tree$maps)
  cumsteps<-c(0,cumsum(nsteps_edges))
  steps<-c()
  names<-c()
  for (i in 1:length(tree$edge[,1])){
    steps<-c(steps,sum(nsteps_edges[1:i])+i-1,sum(nsteps_edges[1:i])+i)
    names<-c(names,(cumsteps[i]+1):cumsteps[i+1],"RM")
  }
  steps<-steps[-length(steps)]
  names<-names[-length(names)]
  custom_cols<-ULT::scale.palette(ncols=max(steps),cols=ord_cols,middle.col=NA,middle=NA,span=c(1,max(steps)),steps=steps[-length(steps)]+1)
  names(custom_cols)<-names
  custom_cols<-custom_cols[-which(names(custom_cols)=="RM")]
  temp_CM$cols<-custom_cols

  for (i in 1:length(tree$edge[,1])){
    names(temp_CM$tree$maps[[i]])<-(cumsteps[i]+1):cumsteps[i+1]
  }

  if(missing(plot.opt)){plot.opt<-NULL}
  do.call(plot,c(list(temp_CM,legend = FALSE),plot.opt))
}
