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
#' @param ncols The number of colors to be outputted. Must be equal or superior to the length of \code{cols}
#' @param cols The colors to be used.
#' @param freqs Optional. The proportions for the repartition of the colors. By default, proportions are all equal.
#' @param steps Optional. The value of the steps the colors represent.
#'
#' @importFrom scales rescale
#' @importFrom stats approx
#'
#' @usage discrete.palette (ncols, cols, freqs, steps)
#'
#' @return
#' A vector that contains \code{ncols} color codes, with all \code{cols} repeated.
#'
#' @export discrete.palette

discrete.palette<-function(ncols,cols,freqs,steps){
  if((missing(freqs)||is.null(freqs)||all(is.na(freqs)))&(missing(steps)||is.null(steps)||all(is.na(steps)))){
    freqs<-rep(ncols%/%length(cols),length(cols))
  }

  if(!missing(steps)){
    discpal<-c()
    if(length(cols)!=length(steps)){
      if(length(steps)<length(cols)){
        warning("Number of steps inferior to the number of colors; rescaled to be all equal")
        steps<-seq(0,1,length(cols))
      }
      if(length(steps)>length(cols)){
        warning("Number of steps superior to the number of colors; only kept the steps corresponding to the given colors")
        steps<-steps[1:length(cols)]
      }
    }
    if(!all(sort(steps)==steps)){
      warning("Steps unordered; automatically reordered them")
      steps<-sort(steps)
    }
    steps<-(steps-min(steps))/(max(steps)-min(steps))
    n_steps<-numeric(length=length(steps))
    for (i in 1:length(n_steps)){
      n_steps[i]<-ifelse(i==length(n_steps),max(steps),mean(c(steps[i+1],steps[i])))-ifelse(i==1,min(steps),mean(c(steps[i-1],steps[i])))
    }
    n_steps<-round(n_steps*ncols)
    if(sum(n_steps)!=ncols){
      med_n_steps<-round(median(1:length(n_steps)))
      n_steps[med_n_steps]<-n_steps[med_n_steps]+1*ifelse(sum(n_steps)>ncols,-1,1)
    }

    for (i in 1:length(n_steps)){
      discpal<-c(discpal,rep(cols[i],n_steps[i]))
    }
  }
  if(!missing(freqs)){
    if(length(cols)!=length(freqs)){
      warning("Number of proportions not equal to the number of colors; rescaled to be all equal")
      freqs<-rep(ncols%/%length(cols),length(cols))
    }
    if(sum(freqs)!=1){
      freqs<-freqs/sum(freqs)
    }
    n_freqs<-round(freqs*ncols)
    discpal<-character(length=sum(n_freqs))

    for (i in 1:length(n_freqs)){
      start<-ifelse(i==1,1,end+1)
      end<-start+n_freqs[i]-1
      discpal[start:end]<-cols[i]
    }
    if(length(discpal)!=ncols){
      discpal<-discpal[round(approx(c(1:length(discpal)),n=ncols)$y-0.499999999999999)]
    }
  }

  return(discpal)
}
