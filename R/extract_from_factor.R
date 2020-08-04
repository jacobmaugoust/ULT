#' @title
#' Extract a factor vector with some levels of an initial factor vector
#'
#' @description
#' This function permits to easily have a reduced factor vector containing some of the original levels
#'
#' @param f1 The initial factor vector
#' @param f2 The vector containing the levels to extract
#'
#' @return
#' A factor vector with only the specified levels in f2
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

  result<-c()
  for (i in 1:length(factor_vector)){
    for (j in 1:length(factor_levels)){
      if(output=="factor"){
        if(factor_vector[i]==factor_levels[j]){
          result<-c(result,as.character(factor_vector[i]))
        }
      }
      if(output=="data"){
        if(factor_vector[i]==factor_levels[j]){
          result<-c(result,x[i])
        }
      }
    }
    if(output=="logical"){
      if(any(factor_vector[i]==factor_levels)){
        result<-c(result,TRUE)
      }
      else{
        result<-c(result,FALSE)
      }
    }
  }
  if(output=="factor"){
    return(as.factor(result))
  }
  else{
    return(result)
  }
}
