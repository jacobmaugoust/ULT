#' @title Map RRphylo evolutionary rates of univariate data on phylogeny
#'
#' @description This function plots the evolutionary rates of a data vector reconstructed with the [RRphylo::RRphylo] method onto a phylogenetic tree.
#' Choices are given regarding the execution of RRphylo, of the colors to be mapped, and to the mapping itself.
#'
#' @details. The goal of this function is to perform separately the phylogenetic comparative method, the color mapping, and the mapping itself, letting the user give additional arguments for each step.
#' The function first runs the [RRphylo::RRphylo] function if a \code{RR} object is not provided, taking into account potential supplementary arguments.
#' Then, the function computes a "basic" mapping using the [phytools::contMap] function but does not show it, taking into account potential additional arguments to be passed to this specific function.
#' Then, the function recreates a color gradient using the [ULT::scale.palette] function, taking into account potential arguments to be passed to this function.
#' Finally, the function replaces the original color gradient by the new, user-specified, one, and plots the new mapping with it. \cr
#' A supplementary interesting point of this function is that the user can provide a covariate and/or a predictor of the variable with only tip values; the function then runs the whole process.
#' Additionally, the user can choose whether to plot the original evolutionary rates or their absolute values.
#'
#' @param RR Optional. The RRphylo object to work with. If not provided, it is computed (which can take time depending on the data).
#' @param tree A phylogenetic tree.
#' @param y A data vector.
#' @param partial.cov Optional. A covariate (see [RRphylo::RRphylo] help page for more details) with data only for tips; full covariate (i.ee., with data for tips and for nodes) is then computed within this function.
#' @param partial.x1 Optional. A predictor (see [RRphylo::RRphylo] help page for more details) with data only for tips; full predictor (i.ee., with data for tips and for nodes) is then computed within this function.
#' @param absrates Logical. Specify whether the user wants to plot natural evolutionary rates (\code{absrates=FALSE}) or absolute evolutionary rates (\code{absrates=TRUE}). Default turned to \code{FALSE}.
#' @param RR.args A list of arguments to be passed to the [RRphylo::RRphylo] function (see function help page for details), with the exception of \code{tree} and \code{y}.
#' @param scale.palette.args A list of arguments to be passed to the [ULT::scale.palette] function (see function help page for details). By default, color is set to a blue-white-red gradient with the white a central color and 0 a central value (\code{scale.palette.args=list(cols=c("blue","white","red"),middle.col="white",middle=0,steps=NA)}) for natural evolutionary rates, and to a heatmap using [RColorBrewer::brewer.pal] (scale.palette.args=list(cols=brewer.pal(9,"YlOrRd"),middle.col=brewer.pal(9,"YlOrRd")[5],middle=NA,steps=NA)) for absolute evolutionary rates.
#' @param contMap.args A list of arguments to be passed to the [phytools::contMap] function (see function help page for details), with the exception of \code{tree}, \code{x}, \code{method}, and \code{anc.states}.
#'
#' @import ape
#' @importFrom RRphylo RRphylo
#' @importFrom phytools contMap
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' require(ape)
#' require(phytools)
#' require(RRphylo)
#' set.seed(1)
#' tree<-rtree(20)
#' y<-fastBM(tree)
#' RR<-RRphylo(tree=tree,y=y)
#' plot.rates.RR(tree=tree,y=y)
#' set.seed(2)
#' x1<-fastBM(tree)*runif(1,-5,5)+runif(1,-5,5)+rnorm(20)
#' plot.rates.RR(tree=tree,y=y,partial.x1=x1)
#' plot.rates.RR(tree=tree,y=y,partial.x1=x1,absrates=TRUE)
#' set.seed(3)
#' cov<-fastBM(tree)*runif(1,-5,5)+runif(1,-5,5)+rnorm(20)
#' plot.rates.RR(tree=tree,y=y,partial.cov=cov)
#' plot.rates.RR(tree=tree,y=y,partial.cov=cov,absrates=TRUE)
#' plot.rates.RR(tree=tree,y=y,partial.x1=x1,partial.cov=cov)
#' plot.rates.RR(tree=tree,y=y,partial.x1=x1,partial.cov=cov,absrates=TRUE)
#'
#' @export plot.rates.RR

plot.rates.RR<-function(RR,tree,y,partial.cov,partial.x1,absrates=FALSE,
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

  rates<-RR$rates[,1]
  if(is.null(scale.palette.args)){
    if(absrates){
      rates<-abs(rates)
      scale.palette.args<-list(cols=brewer.pal(9,"YlOrRd"),middle.col=brewer.pal(9,"YlOrRd")[5],middle=NA,steps=NA)
    }
    else{
      scale.palette.args<-list(cols=c("blue","white","red"),middle.col="white",middle=0,steps=NA)
    }
  }
  rates_tips<-rates[which(names(rates)%in%tree$tip.label)]
  rates_nodes<-rates[which(!names(rates)%in%tree$tip.label)]
  temp_CM<-do.call("contMap",c(list(tree=tree,x=rates_tips,method="user",anc.states=unname(rates_nodes)),contMap.args,plot=FALSE))
  temp_CM$cols<-replace(temp_CM$cols,names(temp_CM$cols),do.call("scale.palette",c(list(ncols=length(temp_CM$cols),span=temp_CM$lims),scale.palette.args)))
  plot(temp_CM)
}
