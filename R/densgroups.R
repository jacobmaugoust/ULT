#' @title Multi-subgroup density lines
#'
#' @description This function plots density lines of subgroups for a given variable
#'
#' @param x The initial variable
#' @param g The groups variable, must be of same length than \code{x}
#' @param cols The colors of the different variable lines. By default using the first colors of the [ULT:contrasting.palette] function
#' @param g.levels.order Optional. The order of the levels for the variable \code{g} if different from the default one
#' @param plot.opt Options (as a list) to be passed to the initial \code{plot} function, that sets the global plot frame. By default set to \code{list(xlab=NA,ylab=NA)}
#' @param curve.opt Options (as a list) to be passed to the \code{lines} function used to draw the density lines. By default set to \code{list(lwd=2)}
#' @param legend Logical indicating whether a legend for each subgroup has to be drawn
#' @param legend.opt Options (as a list) to be passed to the \code{legennd} function if a legend is called. By default set to \code{list("x"="topleft",lwd=2,bty="n",legend=levels(g))}
#'
#' @examples
#' set.seed(1)
#' means<-runif(5,-20,20)
#' sds<-runif(5,0.5,3)
#' x<-c()
#' for(i in 1:5){x<-c(x,rnorm(100,means[i],sds[i]))}
#' g<-as.factor(rep(letters[1:5],each=100))
#' densgroups(x,g) # Default plot
#' densgroups(x,g,legend.opt=list(x="topright")) # Default plot changing location of legend
#' densgroups(x,g, # With some customization
#'            cols=c("red","blue","brown","pink","cyan"),
#'            plot.opt=list(xlab="value",ylab=NA,yaxt="n"),
#'            curve.opt=list(lwd=3,lty=2),
#'            legend.opt=list("x"="topright",lwd=2,lty=2,bty="n",legend=paste0(levels(g),", mu=",round(means,2),", sd=",round(sds,2))))
#'
#' @export

densgroups<-function(x,g,cols,g.levels.order,plot.opt,density.opt,curve.opt,legend=TRUE,legend.opt){
  if(!is.factor(g)){
    if(missing(g.levels.order)){
      g<-as.factor(g)
    }
    else{
      g<-factor(g,g.levels.order)
    }
  }

  max.d<-c()
  for(i in 1:nlevels(g)){
    max.d<-c(max.d,max(hist(x[g==levels(g)[i]],plot=FALSE)$density))
  }

  if(missing(cols)){
    cols<-contrasting.palette(nlevels(g))
  }
  else if(length(cols)<nlevels(g)){
    cols<-c(cols,contrasting.palette(nlevels(g)-length(cols)))
  }

  if(missing(plot.opt)){
    plot.opt<-list(xlab=NA,ylab=NA)
  }
  else{
    plot.opt<-c(plot.opt,
                if(!"xlab"%in%names(plot.opt)){list(xlab=NA)},
                if(!"ylab"%in%names(plot.opt)){list(ylab=NA)})
  }

  if(missing(density.opt)){density.opt<-NULL}

  if(missing(curve.opt)){
    curve.opt<-list(lwd=2)
  }
  else{
    curve.opt<-c(curve.opt,
                 if(!"lwd"%in%names(curve.opt)){lwd=2})
  }

  do.call("plot",c(list(x=1,y=1,type="n",xlim=range(x),ylim=c(0,max(max.d))),plot.opt))

  for(i in 1:nlevels(g)){
    do.call("lines",c(list(x=do.call("density",c(list(x=x[g==levels(g)[i]]),density.opt)),col=cols[i]),curve.opt))
  }

  if(legend){
    if(missing(legend.opt)){
      legend.opt<-list("x"="topleft",lwd=2,bty="n",legend=levels(g))
    }
    else{
      legend.opt<-c(legend.opt,
                    if(!"x"%in%names(legend.opt)){list(x="topleft")},
                    if(!"lwd"%in%names(legend.opt)){list(lwd=2)},
                    if(!"bty"%in%names(legend.opt)){list(bty="n")},
                    if(!"legend"%in%names(legend.opt)){list(legend=levels(g))})
    }
    do.call("legend",c(list(col=cols),legend.opt))
  }
}
