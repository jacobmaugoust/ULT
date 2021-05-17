#' @title
#' Function to provide a user-defined discrete color palette
#'
#' @description
#' This function aims to give a 'discrete' palette, that is, a palette without gradient and with some user-defined repeated colors.
#'
#' @details
#' The function takes as input the colors the user wants to repeat, and optionally their proportions if it is not intended to distribute them equally.
#' The equal distribution may not be fully true: if needed, proportions are rounded to have an integer number of repeats of each color. The goal is rather to output the specified number of colors the user wants.
#'
#' @param ncols The number of colors to be outputted.
#' @param cols The colors to be used.
#' @param prop Optional. The proportions for the repartition of the colors. By default, proportions are all equal.
#'
#' @importFrom scales rescale
#'
#' @usage discrete.palette (ncols, cols, prop)
#'
#' @return
#' A vector that contains \code{ncols} color codes, with all \code{cols} repeated.
#'
#' @export discrete.palette

discrete.palette<-function(ncols,cols,prop){
  d<-ncols%/%length(cols)+1
  r<-abs(ncols-d*length(cols))

  if(missing(prop)){
    prop<-round(rescale(d/100*c(1:length(cols)),c(0,1)),2)
  }
  if(length(cols)!=length(prop)){
    prop<-round(rescale(d/100*c(1:length(cols)),c(0,1)),2)
    warning("Number of proportions not equal to the number of colors; rescaled to be all equal")
  }
  if(max(prop)<1){
    prop<-round(rescale(prop,c(0,1)),2)
    warning("Maximal proportion inferior to one; rescaled to be all equal")
  }

  discpal<-character(length=ncols)

  for (i in 1:(length(prop)-1)){
    start<-ifelse(i==1,1,prop[i]*100+1)
    end<-prop[i+1]*ncols
    discpal[start:end]<-cols[i]
  }

  return(discpal)
}
