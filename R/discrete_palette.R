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
#' @param freqs Optional. The proportions for the repartition of the colors. By default, proportions are all equal.
#' @param steps Optional. The value of the steps the colors represent.
#'
#' @importFrom scales rescale
#'
#' @usage discrete.palette (ncols, cols, freqs, steps)
#'
#' @return
#' A vector that contains \code{ncols} color codes, with all \code{cols} repeated.
#'
#' @export discrete.palette

discrete.palette<-function(ncols,cols,freqs,steps){
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
    n_steps<-round(steps*ncols)
    n_steps[2:length(n_steps)]<-n_steps[2:length(n_steps)]-n_steps[1:(length(n_steps)-1)]
    if(sum(n_steps)!=ncols){
      med_n_steps<-round(median(1:length(n_steps)))
      n_steps[med_n_steps]<-n_steps[med_n_steps]+1*ifelse(sum(n_steps)>ncols,-1,1)
    }

    for (i in 1:length(n_steps)){
      discpal<-c(discpal,rep(cols[i],n_steps[i]))
    }
  }
  if(!missing(freqs)){
    discpal<-character(length=ncols)

    if(missing(freqs)){
      freqs<-rep(ncols%/%length(cols),length(cols))
    }
    if(length(cols)!=length(freqs)){
      warning("Number of proportions not equal to the number of colors; rescaled to be all equal")
      freqs<-rep(ncols%/%length(cols),length(cols))
    }
    if(sum(freqs)!=1){
      freqs<-freqs/sum(freqs)
    }
    n_freqs<-round(freqs*ncols)
    if(sum(n_freqs)!=ncols){
      if(length(which(freqs==min(freqs)))==1){
        med_n_freqs<-which(freqs==min(freqs))
      }
      else{
        med_n_freqs<-round(median(1:length(n_freqs)))
      }
      n_freqs[med_n_freqs]<-n_freqs[med_n_freqs]+1*ifelse(sum(n_freqs)>ncols,-1,1)
    }

    for (i in 1:length(n_freqs)){
      start<-ifelse(i==1,1,end+1)
      end<-start+n_freqs[i]-1
      discpal[start:end]<-cols[i]
    }
  }

  return(discpal)
}
