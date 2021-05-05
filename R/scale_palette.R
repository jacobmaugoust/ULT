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
#' @param cols The colors to be used. At least two colors are expected, three if the user wants to give a 'middle'. Without color given, the assumed vector is c("white","black").
#' @param middle.col The color for the 'center' of the scale. It has to be one of the \code{cols} elements. If not, it is computed as being the middle between the two most extreme given colors, weighted by the \code{middle} of the \code{span} if they are provided.
#' @param span A numeric vector of two values that are the extremities of the range of the values used. For instance, it can be the minimal and maximal values of a continuous numeric vector the colors are intended to follow.
#' @param middle A numeric value that is the 'center' of the color range. This is particularly useful if one wants an asymmetrical range of colors.
#' @param invert Logical. If \code{TRUE}, it flips the color ramp, thus by default set to \code{FALSE}.
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
#' plot(-50:100,-50:100,pch=21,col=NA,bg=scale.palette(ncol=151,cols=c("blue","white","red"),middle.col="white",span=c(-50,100),middle=0))
#'
#' # For a vector following an already defined-gradient, and by skewing its center at 0.75 (i.e., with more reds at the middle):
#' require(RColorBrewer)
#' user_color_gradient<-brewer.pal(8,"YlOrRd") # To have a regular "heat" gradient from yellow to orange to red
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=user_color_gradient,middle.col=user_color_gradient[6],span=c(0,100)))
#' # Or with by skewing the center at 0.25 (i.e., with more yellows):
#' plot(1:100,1:100,pch=21,col=NA,bg=scale.palette(ncols=100,cols=user_color_gradient,middle.col=user_color_gradient[2],span=c(0,100)))
#'
#' @export

scale.palette<-function (ncols,cols,middle.col=NA,span=NA,middle=NA,invert=FALSE){
  if(missing(ncols)){
    stop("No desired number of color provided")
  }
  if(missing(cols)){
    warning("No color provided, palette will be from white to black")
    cols<-c("white","black")
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

  if(is.na(middle)){
    middle<-mean(span)
  }

  old.middle.col<-c()
  if(is.na(middle.col)){
    middle.col<-"white"
    old.middle.col<-NA
  }
  else{
    rm(old.middle.col)
  }

  if(is.na("span")){
    span<-c(0,1)
  }

  middle<-(middle-span[1])/(span[2]-span[1])
  span<-(span-span[1])/(span[2]-span[1])

  if(middle.col%in%cols){
    n_middle_col<-which(middle.col==cols)
  }
  else{
    firstcol<-col2rgb(cols[1])
    lastcol<-col2rgb(cols[length(cols)])
    middle.col<-rgb(red=sqrt(firstcol[1]^2*middle+lastcol[1]^2*(1-middle)),green=sqrt(firstcol[2]^2*middle+lastcol[2]^2*(1-middle)),blue=sqrt(firstcol[3]^2*middle+lastcol[3]^2*(1-middle)),maxColorValue = 255)
    cols<-c(cols[1],middle.col,cols[length(cols)])
    n_middle_col<-which(middle.col==cols)
    if(is.na(old.middle.col)){
      warning("No middle color provided; replacing it by average color between the two extreme color provided")
    }
    else{
      warning("Middle color provided is not in the provided vector color; replacing it by average color between the two extreme color provided")
    }
  }

  final_palette<-c(colorRampPalette(cols[1:n_middle_col])(ncols*(middle-span[1])),colorRampPalette(cols[n_middle_col:length(cols)])(ncols*(span[2]-middle)))

  return(final_palette)
}
