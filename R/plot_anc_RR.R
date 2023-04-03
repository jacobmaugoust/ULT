#' @title Map univariate data on phylogeny following RRphylo method.
#'
#' @description This function plots a data vector onto a phylogenetic tree with ancestral reconstructions following [RRphylo::RRphylo] method.
#' Choices are given regarding the execution of RRphylo, of the colors to be mapped, and to the mapping itself.
#'
#' @details. The goal of this function is to perform separately the phylogenetic comparative method, the color mapping, and the mapping itself, letting the user give additional arguments for each step.
#' The function first runs the [RRphylo::RRphylo] function if a \code{RR} object is not provided, taking into account potential supplementary arguments.
#' Then, the function computes a "basic" mapping using the [phytools::contMap] function but does not show it, taking into account potential additional arguments to be passed to this specific function.
#' Then, the function recreates a color gradient using the [ULT::scale.palette] function, taking into account potential arguments to be passed to this function.
#' Finally, the function replaces the original color gradient by the new, user-specified, one, and plots the new mapping with it. \cr
#' A supplementary interesting point of this function is that the user can provide a covariate and/or a predictor of the variable with only tip values; the function then runs the whole process.
#' Additionally, if a predictor is specified, the user can choose whether to plot the original data or the residuals of the variable vs predictor regression.
#'
#' @param RR Optional. The RRphylo object to work with. If not provided, it is computed (which can take time depending on the data).
#' @param tree A phylogenetic tree.
#' @param y A data vector.
#' @param partial.cov Optional. A covariate (see [RRphylo::RRphylo] help page for more details) with data only for tips; full covariate (i.ee., with data for tips and for nodes) is then computed within this function.
#' @param partial.x1 Optional. A predictor (see [RRphylo::RRphylo] help page for more details) with data only for tips; full predictor (i.ee., with data for tips and for nodes) is then computed within this function.
#' @param res.x1 Logical. If a predictor ('full' or 'partial') is provided, to know whether plotted data should be the \code{y} data vector or the residuals of \code{y} against \code{x1}. Default turned to \code{FALSE}.
#' @param RR.args A list of arguments to be passed to the [RRphylo::RRphylo] function (see function help page for details), with the exception of \code{tree} and \code{y}.
#' @param scale.palette.args A list of arguments to be passed to the [ULT::scale.palette] function (see function help page for details). By default, color is set to a blue-white-red gradient with the white a central color and 0 a central value (\code{scale.palette.args=list(cols=c("blue","white","red"),middle.col="white",middle=0,steps=NA)}).
#' @param contMap.args A list of arguments to be passed to the [phytools::contMap] function (see function help page for details), with the exception of \code{tree}, \code{x}, \code{method}, and \code{anc.states}.
#'
#' @import ape
#' @importFrom RRphylo RRphylo
#' @importFrom phytools contMap
#' @importFrom stats residuals
#'
#' @examples
#' require(ape)
#' require(phytools)
#' require(RRphylo)
#' set.seed(1)
#' tree<-rtree(20)
#' y<-fastBM(tree)
#' RR<-RRphylo(tree=tree,y=y)
#' plot.anc.RR(tree=tree,y=y)
#' set.seed(2)
#' x1<-fastBM(tree)*runif(1,-5,5)+runif(1,-5,5)+rnorm(20)
#' plot.anc.RR(tree=tree,y=y,partial.x1=x1)
#' plot.anc.RR(tree=tree,y=y,partial.x1=x1,res.x1=TRUE)
#'
#' @export plot.anc.RR

plot.anc.RR<-function(RR,tree,y,partial.cov,partial.x1,res.x1=FALSE,
                      RR.args=NULL,
                      scale.palette.args=NULL,
                      contMap.args=NULL){
  if(missing(tree)|missing(y)){
    stop(paste0("Please provide ",ifelse(missing(tree)&missing(y),"a phylogenetic tree and a data vector",
                                         ifelse(missing(tree),"a phylogenetic tree","a data vector"))))
  }
  if(missing(RR)){
    if(missing(RR.args)){
      RR.args<-list(cov=NULL,rootV=NULL,aces=NULL,x1=NULL,aces.x1=NULL,clus=0.5)
    }
    if(!missing(partial.cov)|!missing(partial.x1)){
      if(!missing(partial.cov)){
        RR.args$cov<-c(RRphylo(tree,partial.cov)$aces[,1],partial.cov)
      }
      if(!missing(partial.x1)){
        RR.args$x1<-c(RRphylo(tree,partial.x1)$aces[,1],partial.x1)
      }
    }
    RR<-do.call("RRphylo",c(list(tree=tree,y=y),RR.args))
  }
  if(res.x1){
    all_y<-c(RRphylo(tree,y)$aces[,1],y)
    res_y<-residuals(lm(all_y~RR.args$x1))
    y_tips<-res_y[which(names(res_y)%in%tree$tip.label)]
    y_nodes<-res_y[which(!names(res_y)%in%tree$tip.label)]
  }
  else{
    y_tips<-y
    y_nodes<-RR$aces[,1]
  }
  if(missing(scale.palette.args)){
    scale.palette.args<-list(cols=c("blue","white","red"),middle.col="white",middle=0,steps=NA)
  }

  temp_CM<-do.call("contMap",c(list(tree=tree,x=y_tips,method="user",anc.states=y_nodes),contMap.args,plot=FALSE))
  temp_CM$cols<-replace(temp_CM$cols,names(temp_CM$cols),do.call("scale.palette",c(list(ncols=length(temp_CM$cols),span=temp_CM$lims),scale.palette.args)))
  plot(temp_CM)
}
