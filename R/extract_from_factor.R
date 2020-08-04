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
#' f1<-as.factor(c(rep("A",10),rep("B",5),rep("C",8),rep("A",4)))
#' f2<-c("A","B")
#' extract.from.factor(f1,f2)
#'
#' @export

extract.from.factor<-function(f1,f2){
  if(is.factor(f1)==FALSE){
    f1<-as.factor(f1)
  }
  if(is.character(f2)==FALSE){
    f2<-as.character(f2)
  }
  result<-c()
  for (i in 1:length(f1)){
    for (j in 1:length(f2)){
      if(f1[i]==f2[j]){
        result<-c(result,as.character(f1[i]))
      }
    }
  }
  return(as.factor(result))
}
