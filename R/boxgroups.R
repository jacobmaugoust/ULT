#' @title Multivariate and multi-subgroup boxplots
#'
#' @description This function plots boxplots for various variables associated with ranges for various subgroups found into each variable
#'
#' @details This function first plot the boxplots of the different variables inputted. Then, bars of each of the defined subgroups in each variable are plotted alongside each boxplot, showing desired values describing each subgroup variation within each variable (minimal and maximal values by default, but there can also be other ones detailed below).
#' Graphical options to the \code{boxplot} and \code{points} formulas can be passed to customize the plot.
#' There is also a way to rescale all variables in order to be able to compare the variables and the distribution of each subgroup.
#' \cr
#' There are various ways to describe the variation of each subgroup alongside the global variable boxplot, which are not exclusive. One needs to specify them in the \code{group.plot} argument with a character vector which can contain:
#' \describe{
#'  \item{\code{"min"} and \code{"max"}}{the extreme range, encompassing all points}
#'  \item{\code{"mean"}}{the average value of the whole subgroup}
#'  \item{\code{"median"}}{the median of the whole subgroup}
#'  \item{\code{"sd-X"}}{to obtain points from both sides of the mean by a proportion \code{X} (i.e., a number) of the standard deviation of the whole subgroup}
#'  \item{\code{"q-X"}}{to obtain the \code{X}th percent quantile (i.e., a percentage)}
#'  \item{\code{"all"}}{to show all points of the subgroup}
#' }
#'
#' @param data A matrix or data frame to input.
#' @param groups A vector of same length than the number of rows of \code{data} to distinguish subgroups for each variable.
#' @param box.opt Optional. Plotting options to be passed for the \code{\link[graphics]{boxplot}} function. Must be a list, the names being those of the parameter of this function.
#' @param points.opt Optional. Plotting options to be passed for the \code{\link[graphics]{points}} function. Must be a list, the names being those of the parameter of this function.
#' @param box.width The relative width of the boxes, which also governs the gap between the box and the min-max lines of the variable.
#' @param group.plot The type of points to be drawn to represent the variation of each subgroup in each variable. By default set to the minimal and maximal range of the subgroup to show its entire range.
#' @param group.plot.type The type of plot to be drawn for the plot of each group alongside each boxplot. See \code{plot} argument of the \code{\link[graphics]{plot.default}} function.
#' @param names The names to be displayed under each boxplot + lines. Can be the names of the columns of \code{data}. If nothing is provided, will be alphabetical letters.
#' @param x.gap Gap between the plotted data and the left and right axes. Default set to \code{c(0,0)}; if a single number is provided, automatically meant as left gap only.
#' @param ticks Logical. Should the ticks be displayed above each name? By default, \code{ticks=FALSE}.
#' @param prop Logical. Do the boxplots need to be proportional with each other? By default, natural values are plotted (\code{prop=FALSE}).
#' @param range If \code{prop=TRUE}, the range of values for rescaling the variables. By default set to 0-1.
#'
#' @importFrom scales rescale
#' @importFrom stats quantile
#'
#' @examples
#' A<-runif(30,1,10)
#' B<-rexp(30)
#' C<-rnorm(30,0,1)
#' data<-data.frame(A=A,B=B,C=C)
#' groups<-sample(letters[1:4],30,replace=TRUE)
#' boxgroups(data=data,groups=groups,group.plot = c("all"))
#' boxgroups(data=data,groups=groups,points.opt=list("pch"=21,"col"=c("blue","red","green3","orange"),"bg"=c("blue","red","green3","orange"),"lwd"=2))
#'
#' @export

boxgroups<-function(data,groups,box.opt=NULL,points.opt=NULL,box.width=0.3,group.plot=c("min","max"),group.plot.type=c("b","l"),names=NA,ticks=FALSE,x.gap=c(0,0),prop=FALSE,range=NA){
  groups<-factor(groups)
  ng<-nlevels(groups)
  nvar<-dim(data)[2]
  if(prop){
    if(all(is.na(range))){
      range<-c(0,1)
    }
    for (i in 1:nvar){
      data[,i]<-rescale(data[,i],range,c(min(data[,i],na.rm=TRUE),max(data[,i],na.rm=TRUE)))
    }
  }

  boxwex<-rep(box.width,nvar)
  if(length(x.gap)<2){
    x.gap<-c(x.gap,0)
  }
  xlim<-c(1-box.width*3/2-x.gap[1],nvar+box.width/2+x.gap[2])
  do.call("boxplot",c(list("x"=data,"boxwex"=boxwex,"xlim"=xlim,"xaxt"="n"),box.opt))

  for (i in 1:nvar){
    if(!is.null(points.opt)){
      var.opt<-names(points.opt)[mapply(function(x){length(x)},points.opt)==ng]
    }
    else{
      var.opt<-NULL
    }
    for (j in 1:ng){
      temp_data<-na.omit(data[groups==levels(groups)[j],i])
      y<-c()
      funs<-group.plot[which(group.plot%in%c("min","max","mean","median"))]
      if(length(funs)>0){
        for (k in 1:length(funs)){
          y<-c(y,do.call(funs[k],list("x"=temp_data)))
        }
      }
      if(any(group.plot=="all")){
        y<-c(y,temp_data)
      }
      if(any(mapply(function(x){length(x)},strsplit(group.plot,"-"))>1)){
        concerned<-strsplit(group.plot,"-")[mapply(function(x){length(x)},strsplit(group.plot,"-"))>1]
        fun_types<-mapply(function(x){x[1]},concerned)
        fun_numbers<-as.numeric(mapply(function(x){x[2]},concerned))
        quantiles<-fun_numbers[fun_types=="q"]
        sdvar<-fun_numbers[fun_types=="sd"]
        y<-c(y,as.numeric(quantile(temp_data,quantiles/100)))
        for (sd in min(length(sdvar),1):length(sdvar)){
          y<-c(y,mean(temp_data)-sdvar[sd]*sd(temp_data),mean(temp_data)+sdvar[sd]*sd(temp_data))
        }
      }
      x<-rep(i-box.width/ng*j-box.width/2,length(y))
      if(!is.null(points.opt)&length(var.opt)!=0){
        group.opt<-as.list(mapply(function(x){x[j]},points.opt[var.opt]))
      }
      else{
        group.opt<-NULL
      }
      single.opt<-points.opt[!names(points.opt)%in%var.opt]
      for (k in 1:length(group.plot.type)){
        do.call("points",c(list("x"=x,"y"=sort(y),type=group.plot.type[k]),single.opt,group.opt))
      }
    }
  }
  if(all(is.na(names))){
    if(is.null(colnames(data))){
      names<-letters[1:nvar]
    }
    else{
      names<-colnames(data)
    }
  }
  axis(1,at=seq((1-box.width/2),(nvar-box.width/2),1),tick = ticks,labels = names)
}
