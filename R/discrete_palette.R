#' @title
#' Function to provide a user-defined discrete color palette
#'
#' @description
#' This function aims to give a 'discrete' palette, that is, a palette without gradient and with some user-defined repeated colors.
#'
#' @details
#' The function takes as input the colors the user wants to repeat, and optionally their proportions if it is not intended to distribute them equally.
#' The equal distribution may not be fully true: if needed, proportional samples are rounded to have an integer number of repeats of each color. The goal is rather to output the specified number of colors the user wants.
#'
#' @param ncols The number of colors to be outputted.
#' @param cols The colors to be used.
#' @param prop Optional. The proportions for the repartition of the colors. By default, proportions are all equal.
#' @param prop.type Optional. The type of desired proportion. If \code{prop.type="freq"} (the default), one has to provide a vector of same length than the vector \code{cols} with the frequencies of each color. If \code{prop.type="step"}, one has to provide a vector of the length of \code{cols} + 1 with all steps before to after the considered numeric series (i.e., for two colors, the vector is \code{c(x,y,z)}, having the first color between \code{x} and \code{y} and the second between \code{y} and \code{z}).
#'
#' @importFrom scales rescale
#'
#' @usage discrete.palette (ncols, cols, prop)
#'
#' @return
#' A vector that contains \code{ncols} color codes, with all \code{cols} repeated.
#'
#' @export discrete.palette

discrete.palette<-function(ncols,cols,prop,prop.type="freq"){
  if(prop.type=="step"){
    if((length(cols)+1)!=length(prop)){
      if(length(prop)<length(cols)){
        warning("Number of proportions inferior to the number of colors; rescaled to be all equal")
        prop<-c(0,rep(ncols%/%length(cols),length(cols)))
      }
      else{
        if(length(prop)==length(cols)){
          warning("Number of proportions equal to the number of colors; considered them as frequencies")
          prop.type<-"freq"
        }
        else{
          warning("Number of proportions not equal to the number of colors; only kept the steps corresponding to the given colors")
          prop<-prop[1:(length(cols)+1)]
        }
      }
    }
    if(prop.type=="step"){
      prop<-cumsum(prop)
      prop<-(prop-min(prop))/max(prop)
      n_prop<-round((prop[2:length(prop)]-prop[1:(length(prop)-1)])*ncols)
    }
  }
  if(prop.type=="freq"){
    if(missing(prop)){
      prop<-rep(ncols%/%length(cols),length(cols))
    }
    if(length(cols)!=length(prop)){
      warning("Number of proportions not equal to the number of colors; rescaled to be all equal")
      prop<-rep(ncols%/%length(cols),length(cols))
    }
    if(sum(prop)!=1){
      prop<-prop/sum(prop)
    }
    n_prop<-round(prop*ncols)
  }

  if(sum(n_prop)!=ncols){
    if(length(which(prop==min(prop)))==1){
      med_n_prop<-which(prop==min(prop))
    }
    else{
      med_n_prop<-round(median(1:length(n_prop)))
    }
    n_prop[med_n_prop]<-n_prop[med_n_prop]+1*ifelse(sum(n_prop)>ncols,-1,1)
  }

  discpal<-character(length=ncols)

  for (i in 1:length(n_prop)){
    start<-ifelse(i==1,1,end+1)
    end<-start+n_prop[i]-1
    discpal[start:end]<-cols[i]
  }

  return(discpal)
}
