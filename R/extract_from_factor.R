#' @title
#' Extract a vector with some levels of an initial factor vector
#'
#' @description
#' The main goal of this function is to return a subset or a vector able to subset a data vector linked to a factor vector according to some levels of this factor vector.
#'
#' @param factor_vector The initial factor vector
#' @param factor_levels The vector containing the levels to extract
#' @param x Optional. The data (of any kind) whose values are to be extracted from the factor_vector
#' @param output Optional. The type of output of the function. See Value section below.
#'
#' @return
#' If output is "factor", the function will return a shorter factor vector with the levels specified in the factor_levels argument only.
#' If output is "logical", the function will return a logical vector of same length than the factor_vector specifying if each level is equal to the specified factor_levels. This can be especially useful when subsetting a vector according to multiple factor with (optionally) multiple levels for each.
#' If output is "data", the function will return the data in x if the levels of the factor_vector of the same location are aqual to the specified factor_levels.
#'
#' @examples
#' factor_vector<-as.factor(c(rep("A",10),rep("B",5),rep("A",4),rep("C",6)))
#' factor_levels<-c("B","C")
#' x<-c(rep(10,10),rep(5,5),rep(11,4),rep(9,6))
#'
#' extract.from.factor(factor_vector,factor_levels,output="factor")
#' extract.from.factor(factor_vector,factor_levels,output="logical")
#' extract.from.factor(factor_vector,factor_levels,x,output="data")
#'
#' @export

extract.from.factor<-function(factor_vector,factor_levels,x=NA,output){
  if(is.factor(factor_vector)==FALSE){
    factor_vector<-as.factor(factor_vector)
  }
  if(is.character(factor_levels)==FALSE){
    factor_levels<-as.character(factor_levels)
  }
  if(missing(output)){
    if(missing(x)){
      output<-"factor"
    }
    else{
      output<-"data"
    }
  }
  if(length(x)!=length(factor_vector)){
    errorCondition("Error: data vector does not have the same length than the factor vector")
  }
  if(missing(x)==FALSE&output!="data"){
    warning("You provided a data vector but asked for a logical or factor output / did not asked for a data output. It is suggested to add output='data' in the function")
  }

  result<-c()
  for (i in 1:length(factor_vector)){
    if(is.na(factor_vector[i])){
      result<-c(result,FALSE)
    }
    else if(any(factor_vector[i]==factor_levels)){
      result<-c(result,TRUE)
    }
    else{
      result<-c(result,FALSE)
    }
  }

  if(output=="logical"){
    return(result)
  }
  if(output=="factor"){
    return(factor_vector[result])
  }
  if(output=="data"){
    return(x[result])
  }
}
