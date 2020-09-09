#' @title
#' Function to plot two y-variables relative to a x-variable
#'
#' @description
#' This function performs a classic plot of a y-variable (y) relative to a x-variable then superimpose a second plot of a second y-variable (y2) relative to the same x-variable
#'
#' @details
#' The function basically needs three vectors to be performed: the "x", "y" and "y2" variables.
#' It also needs the packages 'scales'.
#' The x-variable can be of any type of vector. Both y and y2 variables must be numeric or integer.
#' This function firstly perform a plot of x and y (by default, "plot" is used with x and y successively, not a ~ formula) then performs a plot of x and of a y-scaled y2 variable (by default, "points" is used).
#' The first y axis is plotted to the left, while the second one is plotted to the right.
#' The second/right y axis is of same size than the first/left y axis, but of different scale.
#' Plot types can be modified, as well as graphical parameters of the plots and of the two y-axis.
#'
#' @param x The x-variable to be plotted. Can be either numeric, integer but also factor or character.
#' @param y The first y-variable to be plotted (y-axis is plotted to the left). Can be either numeric or integer.
#' @param y2 The second y-variable to be plotted (y2-axis is plotted to the right). Can also be either numeric or integer.
#' @param xy.plot.type The plot function used to plot the x-y relationship: can be "plot" or "boxplot". By default, the formula is plot(x,y).
#' @param xy2.plot.type The plot function used to plot the x-y2 relationship: can be "points" or "boxplot". By default, the formula is points(x,y).
#' @param xy.plot.param The graphical parameters to use for the x-y plot; must be a list.
#' @param xy2.plot.param The graphical parameters to use for the x-y2 plot; must be a list.
#' @param y.axis.param The graphical parameters to use for the first/left y axis; must be a list.
#' @param y2.axis.param The graphical parameters to use for the second/right y axis; must be a list.
#'
#' @importFrom scales rescale
#' @importFrom graphics axis
#' @importFrom stats lm
#'
#' @examples
#'
#' # For all three variables being numeric:
#' x<-sort(runif(100,0,1))
#' y<-sort(rnorm(100,0,1))
#' y2<-sort(rexp(100))
#' plot_twosided(x,y,y2,xy.plot.param=list(type="l",col="red",lwd=2),xy2.plot.param=list(type="l",col="blue",lwd=2),y.axis.param=list(col="red"),y2.axis.param=list(col="blue"))
#'
#' # For character (= factorial) x:
#' x<-as.factor(sort(rep(letters[1:10],10)))
#' y<-rexp(100,0.5)
#' y2<-sort(rep(runif(10,15,25),10))
#' plot_twosided(x,y,y2,xy.plot.param=list(col="red",border="blue"),xy2.plot.param=list(type="l",col="green3",lwd=2),y.axis.param=list(col="red"),y2.axis.param=list(col="green3"))
#'
#' @export
plot.twosided<-function(x,y,y2,xy.plot.type,xy2.plot.type,xy.plot.param,xy2.plot.param,y.axis.param,y2.axis.param){
  if(is.character(x)){x<-as.factor(x)}

  if(missing(xy.plot.param)){xy.plot.param<-list(yaxt="n",xlab="",ylab="")}
  else{
    if("xlab"%in%names(xy.plot.param)&"ylab"%in%names(xy.plot.param)){xy.plot.param<-c(yaxt="n",as.list(xy.plot.param))}
    else{
      if(xor("xlab"%in%names(xy.plot.param),"ylab"%in%names(xy.plot.param))){
        if("xlab"%in%names(xy.plot.param)){xy.plot.param<-c(yaxt="n",ylab="",as.list(xy.plot.param))}
        if("ylab"%in%names(xy.plot.param)){xy.plot.param<-c(yaxt="n",xlab="",as.list(xy.plot.param))}
      }
      else{
        xy.plot.param<-c(yaxt="n",xlab="",ylab="",as.list(xy.plot.param))
      }
    }
  }
  if(missing(xy.plot.type)){xy.plot.type<-"plot"}
  do.call(xy.plot.type,c(
    if(xy.plot.type=="boxplot"){list(formula=y~x)}
    else{list(x=x,y=y)},
    xy.plot.param))

  y2_scaled<-rescale(y2,to=c(min(y,na.rm=TRUE),max(y,na.rm=TRUE)),na.rm=TRUE)

  if(missing(xy2.plot.param)){xy2.plot.param<-list(yaxt="n")}
  else{
    if("xlab"%in%names(xy2.plot.param)&"ylab"%in%names(xy2.plot.param)){xy2.plot.param<-c(yaxt="n",as.list(xy2.plot.param))}
    else{
      if(xor("xlab"%in%names(xy2.plot.param),"ylab"%in%names(xy2.plot.param))){
        if("xlab"%in%names(xy2.plot.param)){xy2.plot.param<-c(yaxt="n",ylab="",as.list(xy2.plot.param))}
        if("ylab"%in%names(xy2.plot.param)){xy2.plot.param<-c(yaxt="n",xlab="",as.list(xy2.plot.param))}
      }
      else{
        xy2.plot.param<-c(yaxt="n",xlab="",ylab="",as.list(xy2.plot.param))
      }
    }
  }
  if(missing(xy2.plot.type)){xy2.plot.type<-"points"}
  do.call(xy2.plot.type,c(
    if(xy2.plot.type=="boxplot"){list(formula=y2_scaled~x)}
    else{list(
      if(is.numeric(x)|is.integer(x)){x=x}
      else{x=as.numeric(as.factor(x))},
      y=y2_scaled)},
    if(xy2.plot.type=="boxplot"){c(add=TRUE,xy2.plot.param)}
    else{xy2.plot.param}))

  y_scale<-axis(2,labels=FALSE,tick=FALSE)
  y2_scale<-rescale(y_scale,to=c(min(y2,na.rm=TRUE),max(y2,na.rm=TRUE)))

  y_bar<-pretty(y_scale)
  y_at<-y_scale

  y2_bar<-pretty(y2_scale)
  scaling<-suppressWarnings(summary(lm(y_scale~y2_scale))[[4]])
  y2_at<-y2_bar*scaling[2,1]+scaling[1,1]
  if(y2_at[1]<0){
    y2_bar<-y2_bar[-1]
    y2_at<-y2_at[-1]
  }
  last<-length(y2_bar)
  if(y2_at[last]>max(y_scale)){
    y2_bar<-y2_bar[-last]
    y2_at<-y2_at[-last]
  }

  if(missing(y.axis.param)){
    axis(2,at=y_scale,labels=y_scale)
  }
  else{
    do.call(axis,c(list(side=2,at=y_scale,labels=y_scale),y.axis.param))
  }

  if(missing(y2.axis.param)){
    axis(4,at=y2_at,labels=as.character(y2_bar))
  }
  else{
    do.call(axis,c(list(side=4,at=y2_at,labels=as.character(y2_bar)),y2.axis.param))
  }
}
