#' @title
#' Get color palette depending on frequencies within data of given set of colors
#'
#' @description
#' This function uses a given set of colors and frequencies of them within given dataset and range to derive a set of a given number of colors fitting with the inputted frequencies.
#'
#' @details
#' This function mainly relies on \link[ULT]{discrete.palette} to resize the inputted colors taking into account their frequencies.
#' The frequencies \code{freqs} can be given explicitly (if specifying numeric values), otherwise they are considered to be all equal (if specifying a character).
#' Numeric frequencies can be either proportions or bounds between \code{lims}. In both cases, they depend on a given set of \code{values} and on limits \code{lims} to consider (by default the range of the dataset).
#' If frequencies are proportions, they are converted to bounds to determine the effective width between them (or to one of the limits) relative to the inter-limits width. Hence, bounds can also be directly provided, if known.
#' Such relative widths then serve as the \code{freqs} argument of \link[ULT]{discrete.palette}.
#' If proportions are provided (explicit numeric or just character), it is also possible to say with \code{type} whether these proportions are between-limits width proportion or \code{values} population proportion. In the latter case, quantiles corresponding to the given frequencies \code{freqs} are hence used as bounds.
#' The output will use adjusted frequencies (given the inputted one and the width/proportion choice if applicable) to resize the inputed colors \code{cols} as a new vector of colors of size \code{ncols}
#'
#' @param cols The vector of colors to be resized given data and frequencies
#' @param ncols The number of desired colors for the new vector of colors with adjusted color proportions
#' @param freqs The inputted frequencies of the colors. Can be a character (\code{"equal"} or \code{"even"}, but actually any character, giving equal frequencies to all colors) or a numeric (specifying proportions or bounds to use)
#' @param lims The limits to use to compute bounds (and then adjusted frequencies) for each color. Set by default to the range of \code{values}
#' @param type If specifying proportions to \code{freqs}, whether they represent fractions of the width between \code{lims} (\code{type="width"}) or real proportions of the population of \code{values} (\code{type="proportion"})
#' @param values If \code{type="proportion"}, the data to use to get the real population proportions and then the adjusted color frequencies
#'
#' @importFrom stats approx
#'
#' @export freq.cols

freq.cols<-function(cols,ncols,freqs,lims,type,values){
  quantile.to.freqs<-function(q,lims,length.out){
    if(length(q)==length.out){
      q<-q[-which.max(q)]
    }
    freqs<-(c(q,max(lims))-c(min(lims),q))/diff(lims)
    freqs
  }

  colength<-length(cols)

  if(length(type)>1){
    type<-type[1]
  }

  if(is.character(freqs)&&freqs%in%c("equal","even")){
    freqs<-rep(1/colength,colength)
  }

  if(sum(freqs)!=1){
    freqs<-quantile.to.freqs(freqs,lims,colength)
  }
  else{
    if(length(freqs)!=colength){
      new.freqs<-rep(freqs/round(colength/length(freqs)),each=round(colength/length(freqs)))
      if(length(new.freqs)!=colength){
        diff<-abs(length(new.freqs)-colength)
        while(diff!=0){
          freq.order<-order(freqs,decreasing=TRUE)
          for(i in freq.order){
            i_freq<-which(new.freqs==(freqs[i]/round(colength/length(freqs))))[1]
            if(length(new.freqs)<colength){
              new.freqs<-c(new.freqs[1:i_freq],new.freqs[i_freq],new.freqs[(i_freq+1):length(new.freqs)])
            }
            else{
              new.freqs<-new.freqs[-i_freq]
            }
            diff<-abs(length(new.freqs)-colength)
            if(diff==0){break}
          }
        }
      }
      freqs<-new.freqs
      if(sum(freqs)!=1){
        freqs<-freqs/sum(freqs)
      }
    }
    if(type=="proportion"){
      freqs<-quantile.to.freqs(quantile(values,cumsum(freqs)),lims,colength)
    }
  }

  all.cols<-do.call("discrete.palette",list("cols"=cols,
                                            "ncols"=max(ncols,round(max(1/freqs))),
                                            "freqs"=freqs))

  cols<-all.cols[round(approx(c(1:length(all.cols)),n=ncols)$y-0.499999999999999)]
  return(cols)
}
