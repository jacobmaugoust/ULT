#' @title
#' Binarize a vector
#'
#' @description
#' This function transforms a vector to a vector composed of zeros and ones.
#'
#' @details
#' For given zero and one 'templates', the function transforms each one in zero or in one and return a 0/1 vector.
#'
#' @param x The vector to binarize
#' @param zero The pattern of zeros to be find in x. By default, set to NA. Can be anything.
#' @param one The pattern of ones to be find in \code{x}. By default, set to \code{"else"}, which means all that is not of the \code{zero} pattern. Otherwise, can be anything.
#' @param drop If there are values not fitting both previous patterns, theses values are returned to NA. By default (\code{drop=TRUE}), these NA are left in the output vector. If you want to discard them, please specify \code{drop=TRUE}.
#' @param output The output of the binarization. By default to 0/1, but it can be changed with a vector or a list of length 2.
#'
#' @return
#' The binarized vector of the input vector, with only 0 and 1 (and NA if applicable).
#'
#' @examples
#'
#' set.seed(1)
#' x<-runif(100,0,1) # Create a vector of random values
#' x[as.integer(runif(10,1,100))]<-NA # Randomly replace some of these by NA's
#' binarize(x) # Binarize it with default option
#' all(which(is.na(x))==which(binarize(x)==0)) # Check if all NA are zero
#' all(which(!is.na(x))==which(binarize(x)==1)) # Check if all non-NA are one
#'
#' x<-c(rep(0,20),2,rep(1,20)) # Create a vector with three repeated numbers: 0, 1, and 2
#' binarize(x,zero=0,one=1) # Binarize it by saying that zeros are 0, and ones are 1 (thus, two is not allocated, and it is replaced by NA)
#' binarize(x,zero=0,one=1,drop=TRUE) # Dropping the NA value (which are the no 0 or 1)
#' length(x) # How many numbers are in x
#' length(binarize(x,zero=0,one=1)) # How many numbers if we binarize with default (normally, same length than x, because the 2 is replaced by NA)
#' length(binarize(x,zero=0,one=1,drop=FALSE)) # If we decide do don't drop, the NA is still there, but there is no warnings
#' length(binarize(x,zero=0,one=1,drop=TRUE)) # If we decide to drop, the NA is removed, and the length is shorter of one
#'
#' @export

binarize<-function(x,zero=NA,one="else",drop,output=c(0,1)){
  binarized<-rep(NA,length(x))

  if(all(is.null(zero))){
    zeros<-which(is.null(x))
  }
  if(all(is.na(zero))){
    zeros<-which(is.na(x))
  }
  else{
    zeros<-which(x%in%zero)
  }

  if(length(one)==1&&one=="else"){
    ones<-which(x%in%x[-zeros])
  }
  else{
    ones<-which(x%in%one)
  }

  if(missing(drop)){
    drop<-FALSE
    if(any((x%in%x[zeros]|x%in%x[ones])==FALSE)){
      warning("There are values in x that not match the zero nor the one patterns; they have been treated as NA and left in the resulting vector. Please specify drop=TRUE in the function if you want to discard them.")
    }
  }

  if(length(output)>2){
    warning("There are more than two values in the possible output; please either provide two values, or a list of two vectors; yet, only the two first values provided are used")
    output<-output[1:2]
  }

  if(!is.list(output)){
    output<-list(output[1],output[2])
  }

  if(length(output[[1]])==1){
    binarized[zeros]<-output[[1]]
  }
  else{
    if(length(zeros)!=length(output[[1]])){
      warning("less output values than zero's")
    }
    suppressWarnings(binarized[zeros]<-output[[1]])
  }

  if(length(output[[2]])==1){
    binarized[ones]<-output[[2]]
  }
  else{
    if(length(ones)!=length(output[[2]])){
      warning("less output values than one's")
    }
    suppressWarnings(binarized[ones]<-output[[2]])
  }

  if(drop){
    the_NAs<-which((x%in%x[zeros]|x%in%x[ones])==FALSE)
    binarized<-binarized[-the_NAs]
  }

  return(binarized)
}
