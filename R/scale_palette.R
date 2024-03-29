#' @title
#' Function to provide a user-defined scaled color palette
#'
#' @description
#' This function gives a vector of colors to use following given input colors, with a range following a given span, and with an adjustable 'center' of the scale (i.e., the palette can be asymmetric).
#' This is particularly useful for color gradients that follow continuous, numeric, vectors: the user can thus adjust which color is in the 'center' of the scale, and what are the 'color paths' on each side of that 'center'.
#'
#' @details
#' The function uses the function \code{colorRampPalette} in order to provide a continuous color gradient following a given input of colors and with a given length.
#' This function enables the user to get an asymmetrical palette, either in the color repartition or in the extreme colors of the palette.
#'
#' @param ncols The number of colors to be outputted. It is scaled linearly according to the given \code{span}.
#' @param cols Optional. The colors to be used. At least two colors are expected, three if the user wants to give a 'middle'. Without color given, the assumed vector is \code{c("white","black")}.
#' @param middle.col Optional. The color for the 'center' of the scale. It has to be one of the \code{cols} elements. If not, it is computed as being the middle color between the two most extreme given colors, weighted by the \code{middle} of the \code{span} if they are provided. Type \code{middle.col=NA} if no middle color is desired but is acknowledged, avoiding warning messages.
#' @param span Optional. A numeric vector of two values that are the extremities of the range of the values used. For instance, it can be the minimal and maximal values of a continuous numeric vector the colors are intended to follow. If not provided, it is set to \code{c(0,1)}. Type \code{span=NA} if no span is desired but is acknowledged, avoiding warning messages.
#' @param middle Optional. A numeric value that is the 'center' of the color range. This is particularly useful if one wants an asymmetrical range of colors. If not provided, it is the half of the span. Type \code{middle=NA} if a (default) halfway middle is desired but is acknowledged, avoiding warning messages.
#' @param steps Optional. A vector giving the steps between colors if one wants a multi-step scaled palette. Must be the length of the length of \code{cols} minus two (i.e., the two extreme colors).
#' @param invert Optional. If \code{invert = TRUE}, it flips the color ramp; by default set to \code{FALSE}.
#'
#' @importFrom grDevices colorRampPalette
#'
#' @usage scale.palette(ncols, cols, middle.col, span, middle, steps, invert=FALSE)
#'
#' @return
#' A vector that contains \code{ncols} color codes.
#'
#' @examples
#' # For a vector going from 0 to 100, with colors from red to blue:
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=c("red","blue"),span=c(0,100)))
#'
#' # The same example but with a 'middle color' at the 20, thus skewing the color ramp with more reds between 0 and 20 and more blues between 20 and 100:
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=c("red","blue"),span=c(0,100),middle=20))
#'
#' # For a vector going from 0 to 100, with colors from red to blue to green:
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=c("red","green3","blue"),middle.col="green3",span=c(0,100)))
#'
#' # The same example but with a 'middle color' at the 20, thus skewing the color ramp with reds to blues between 0 and 20 and blues to greens between 20 and 100:
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=c("red","green3","blue"),middle.col="green3",span=c(0,100),middle=20))
#'
#' # For a vector going from -50 to 100, with blues towards -50, reds towards 100, and a white center at 0:
#' plot(-50:100,-50:100,pch=21,col=NA,bg=scale.palette(ncols=151,cols=c("blue","white","red"),middle.col="white",span=c(-50,100),middle=0))
#'
#' # For a vector following an already defined-gradient, and by skewing its center at 0.75 (i.e., with more reds at the middle):
#' require(RColorBrewer)
#' user_color_gradient<-brewer.pal(8,"YlOrRd") # To have a regular "heat" gradient from yellow to orange to red
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=user_color_gradient,middle.col=user_color_gradient[6],span=c(0,100)))
#' # Or with by skewing the center at 0.25 (i.e., with more yellows):
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=user_color_gradient,middle.col=user_color_gradient[2],span=c(0,100)))
#'
#' # For a multi-steps gradient
#'
#' require(ULT)
#' ncols<-round(runif(1,1000,10000))
#' cols<-contrasting.palette(1:round(runif(1,2,18)))
#' span<-c(-13,17.5)
#' steps<-seq(span[1],span[2],length.out=(length(cols)))[-c(1,length(cols))]
#'
#' bg_cols<-scale.palette(ncols=ncols,cols=cols,middle.col=NA,span=span,middle=NA,steps=steps,invert=FALSE)
#' plot(seq(span[1],span[2],length.out=ncols),seq(span[1],span[2],length.out=ncols),pch=21,col=NA,bg=bg_cols)
#' abline(v=span[1],lwd=2,col=cols[1])
#' for (i in 1:length(steps)){
#'   abline(v=steps[i],lwd=2,col=cols[i+1])
#' }
#' abline(v=span[2],lwd=2,col=cols[length(cols)])
#'
#' @export scale.palette

scale.palette<-function(ncols,cols,middle.col,span,middle,steps,invert=FALSE){
  if(missing(ncols)){
    stop("No desired number of color provided")
  }
  if(missing(cols)||is.null(cols)||all(is.na(cols))){
    if(missing(cols)){
      warning("No color provided, palette will be from white to black")
    }
    cols<-c("white","black")
  }
  if(missing(span)||is.null(span)||all(is.na(span))){
    if(missing(span)){
      warning("No range provided to scale the color gradient; 0-1 span taken by default")
    }
    span<-c(0,1)
  }
  if(missing(middle)||is.null(middle)||all(is.na(middle))){
    if(missing(middle)){
      warning("No middle provided to adjust a middle color; half of the span taken by default")
    }
    middle<-mean(span)
  }
  if(missing(middle.col)||is.null(middle.col)||all(is.na(middle.col))){
    if(missing(middle.col)){
      warning(paste0("No middle color provided; replacing it by the ",if(length(cols)==2||length(cols)%%2==0){"average color between the two "},if(length(cols)==2){"provided"}else{"closest"}," color",if(length(cols)%%2==0){"s"},if(length(cols)>2){" to the middle of the colors vector"}))
    }

    if(length(cols)==2){
      middle.col<-rgb(t(apply(col2rgb(cols),1,function(x){sqrt(x[1]^2*middle/diff(span)+x[2]^2*(1-middle/diff(span)))})),maxColorValue = 255)
    }
    else{
      if(((length(cols)-2)%%2)==1){
        middle.col<-cols[(length(cols)+1)/2]
      }
      else{
        middle.col<-rgb(t(apply(col2rgb(cols[length(cols)/2+c(0,1)]),1,function(x){sqrt(x[1]^2*middle/diff(span)+x[2]^2*(1-middle/diff(span)))})),maxColorValue = 255)
      }
    }
    if(!middle.col%in%cols){
      cols<-c(cols[1:(length(cols)/2)],middle.col,cols[length(cols)/2+(1:(length(cols)/2))])
    }
  }
  if(missing(steps)||is.null(steps)||all(is.na(steps))||length(steps)!=(length(cols)-2)){
    if(!(missing(steps))&&(is.null(steps)||all(is.na(steps))&&length(steps)==(length(cols)-2))){
      warning("Provided steps are not equal to the number of non-extreme colors; converted steps to NA")
    }
    steps<-numeric(length(cols)-2)
    n.middle.col<-which(cols==middle.col)-1
    steps[1:n.middle.col]<-seq(span[1],middle,length.out=length(c(1:n.middle.col))+1)[-1]
    steps[n.middle.col:length(steps)]<-seq(middle,span[2],length.out=length(c(n.middle.col:length(steps)))+1)[-(length(c(n.middle.col:length(steps)))+1)]
  }

  if(is.matrix(cols)==TRUE){
    if(dim(as.matrix(cols))[1]==3){
      if((names(cols[1,])=="red"&names(cols[2,])=="green"&names(cols[3,])=="blue")||all(rownames(cols)==c("red","green","blue"))){
        if(any(cols>1)){
          cols<-cols/255
        }
        else{
          if(all(cols>=0)){
            cols<-cols
          }
        }
        old_cols<-cols
        cols<-c(NA)
        for (i in 1:dim(as.matrix(old_cols))[2]){
          cols[i]<-rgb(red=old_cols[1,i],green=old_cols[2,i],blue=old_cols[3,i])
        }
      }
      else{
        warning("Front color is not a RGB matrix, the red/green/blue row names are missing as the output of the col2rgb function; try something like col2rgb(rgb(t(front_color))) or name the rows")
      }
    }
    else{
      warning("Front color has not the size of a RGB matrix as the output of the col2rgb function")
    }
  }
  else{
    if(is.character(cols)==TRUE){
      cols<-cols
    }
    else{
      warning("Front color is not a RGB code color, a HEX code color or a named color")
    }
  }

  if(invert){
    cols<-rev(cols)
  }

  final_palette<-c()

  if(diff(span)==0){
    final_palette<-rep(middle.col,ncols)
  }
  else{
    steps<-c(span[1],steps,span[2])
    steps<-(steps-steps[1])/diff(span)
    n_steps<-round(ncols*steps)
    for (i in 1:(length(steps)-1)){
      if(i==1){local_palette<-c(colorRampPalette(cols[i:(i+1)])(n_steps[i+1]-n_steps[i]))}
      else{local_palette<-c(colorRampPalette(cols[i:(i+1)])(n_steps[i+1]-n_steps[i]+1))[-1]}
      final_palette<-c(final_palette,local_palette)
    }
  }

  return(final_palette)
}
