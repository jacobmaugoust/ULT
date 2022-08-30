#' @title Draw arrows next to phylogenetic taxa
#'
#' @description By giving taxon numbers of a phylogenetic tree and graphical parameters, the function plots the phylogeny (if not already plotted) and arrows next to taxa for graphical purposes.
#'
#' @param phy A phylogenetic tree.
#' @param tax A vector with the name or the number of the taxa (tips and/or nodes) for plotting the arrows next to.
#' @param side A vector to specify the side for plotting the arrows. Can be one \code{1} or \code{"topleft"} (the default), \code{2} or \code{"topright"}, \code{3} or \code{"bottomleft"}, or \code{4} or \code{"bottomright"}.
#' @param angle The angle of the arrow in degrees. Can be a value between 0 and 360, somewhat overriding the \code{side} parameter. Default set to 45.
#' @param length The length of the arrow. Default is to one tenth of the x span, but the user may try several values since this length if then expressed using both x and y.
#' @param offset The length of the distance between the taxon and the end of the arrow. Default is to one tenth of the length, but the user may try several values since this length if then expressed using both x and y.
#' @param arrow.fun The function to be used to draw the arrows, keeping in mind that the system used is here by providing x0/x1/y0/y1 coordinates. Default to \code{arrows}, but the user may want to use the [shape::Arrows] function.
#' @param arrow.pars A list of arguments to be passed to the arrow function used. See [base::arrows] for the base function arguments; see [shape::Arrows] for the function from the \code{shape} package.
#'
#' @import ape
#' @importFrom phytools nodeHeights
#'
#' @examples
#' require(ape)
#' set.seed(1)
#' phy<-rtree(20)
#' plot(phy)
#' tax<-c(2,7,32)
#' tax.arrows(phy,tax)
#' tax.arrows(phy,tax,angle=0)
#' tax.arrows(phy,tax,angle=90)
#' tax.arrows(phy,tax,angle=45,arrow.pars=list(length=0.1,angle=30))
#' require(shape)
#' tax.arrows(phy,tax,angle=45,arrow.pars=list(arr.length=0.5,arr.type="triangle"),arrow.fun="Arrows")
#' tax.arrows(phy,tax,angle=45,arrow.pars=list(arr.length=0.2,arr.type="triangle",arr.col="red"),arrow.fun="Arrows")
#'
#' @export tax.arrows


tax.arrows<-function(phy,tax,side="topleft",angle=45,length,offset,arrow.fun="arrows",arrow.pars){
  NH<-nodeHeights(phy)
  x<-NH[,2]
  x.coo<-c(x[phy$edge[,2]<=Ntip(phy)],0,x[phy$edge[,2]>Ntip(phy)])
  y.coo<-node.height(phy)

  if(is.character(side)){
    side<-which(side==c("topleft","topright","bottomleft","bottomright"))
  }
  else if(side>4){
    stop("You did not provided a valid side; please provide a value between 1 and 4 or a character (topleft, topright, bottomleft, bottomright)")
  }

  suppressWarnings(try(par(new=TRUE)))
  check.plot<-suppressWarnings(try(par(new=TRUE)))$new
  if(!check.plot){
    plot(phy)
  }

  xy_ratio<-diff(par("usr")[3:4])/diff(par("usr")[1:2])*diff(grconvertX(par("usr")[1:2],"user","inches"))/diff(grconvertY(par("usr")[3:4],"user","inches"))

  if(missing(length)){
    length<-max(NH)/10
  }
  x.L<-length*cos(angle*pi/180)
  y.L<-length*sin(angle*pi/180)*xy_ratio

  if(missing(offset)){
    offset<-length/10
  }
  x.off<-offset*cos(angle*pi/180)
  y.off<-offset*sin(angle*pi/180)*xy_ratio

  if(missing(arrow.pars)){
    arrow.pars<-NULL
  }
  for(i in 1:length(tax)){
    if(is.character(tax[i])){
      who<-which(tax[i]==c(phy$tip.label,phy$node.label))
    }
    else{
      who<-tax[i]
    }
    x0<-x.coo[who]+ifelse(side==1|side==3,-1,1)*x.off+ifelse(side==1|side==3,-1,1)*x.L
    x1<-x.coo[who]+ifelse(side==1|side==3,-1,1)*x.off
    y0<-y.coo[who]+ifelse(side==1|side==3,1,-1)*y.off+ifelse(side==1|side==3,1,-1)*y.L
    y1<-y.coo[who]+ifelse(side==1|side==3,1,-1)*y.off
    do.call(arrow.fun,c(list(x0=x0,x1=x1,y0=y0,y1=y1),arrow.pars))
  }
}
