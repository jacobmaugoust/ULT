#' @title T-test of two regression slopes or intercepts
#'
#' @description
#' This function performs a t-test to test whether two regression slopes or two regression intercepts are equal or not.
#'
#' @details
#' The t-test can be performed on two x and y value vectors if also bringing a factor vector separating the x and y data in two subgroups.
#' It can also be performed on two independant x and y value vectors, that are x, y, x2 and y2.
#' Finally, it can be performed on two lm outputs (no matter if the output comes from \code{lm()} or \code{summary(lm()))}.
#'
#' @param x the x-values to be used. Can contain the two groups of data if the parameter \code{factor} is provided. Otherwise, contain the x-values for the first group.
#' @param y the y-values to be used. As for \code{x}, it can contain the two groups of data if the parameter \code{factor} is provided. Otherwise, contain the y-values for the first group.
#' @param factor mandatory if providing \code{x} and \code{y} data to separate two groups.
#' @param x2 the x-values to be used for the second group.
#' @param y2 the y-values to be used for the second group.
#' @param object_1 the output of the linear model for the first group; can be an output of \code{"lm"} or of \code{"summary(lm)"}
#' @param object_2 the output of the linear model for the second group; can be an output of \code{"lm"} or of \code{"summary(lm)"}
#' @param parameter the parameter of the regression line to be tested; by default the \code{"slope"}, otherwise can be the \code{"intercept"}
#'
#' @return Returns an object containing the absolute difference between the two estimated parameters ($Estimate), the standard error of this difference ($Std. Error), the t-value of the t-test ($t value), the degrees of freedom of the test ($df) and the p-value of the bilateral t-test ($p-value).
#'
#' @examples
#'
#' # For an expected equality of slopes :
#' x_A<-sort(runif(50,0,10))
#' x_B<-c()
#' for (i in 1:50){
#'   x_B[i]<-x_A[i]+runif(1,-0.1,0.1)
#' }
#' x<-c(x_A,x_B)
#' # Otherwise :
#' x<-c(sort(runif(50,0,10)),sort(runif(50,0,10)))
#'
#' y<-c(sort(runif(50,0,2)),sort(runif(50,2,4)))
#' fac<-as.factor(c(rep("A",50),rep("B",50)))
#' # Using the x, y and fac vectors only :
#' t.slope.test(x=x,y=y,factor=fac)
#' # Using the x, y, x2 and y2 parameters :
#' t.slope.test(x=x[fac=="A"],y=y[fac=="A"],x2=x[fac=="B"],y2=y[fac=="B"])
#' # Using the output_1 and output_2 parameters
#' lm_1<-lm(y[fac=="A"]~x[fac=="A"])
#' lm_2<-lm(y[fac=="B"]~x[fac=="B"])
#' t.slope.test(output_1=lm_1,output_2=lm_2)
#'
#' @export

t.slope.test<-function(x,y,factor,x2,y2,object_1,object_2,parameter="slope"){
  if((missing(object_1)==FALSE|missing(object_2)==FALSE)&(missing(x)==FALSE|missing(y)==FALSE|missing(factor)==FALSE|missing(x2)==FALSE|missing(y2)==FALSE)==FALSE){
    errorCondition("You provided a lm output and a vector; please choose the method to use and provide either two lm outputs or x and y data together with a factorial vector distinguishing two groups")
  }
  if(missing(object_1)&missing(object_2)){
    if(missing(x2)&missing(y2)){
      if(missing(x)){
        if(missing(y)){
          errorCondition('Only the vector "factor" is provided; please also provide "x" and "y" vectors or provide two lm outputs')
        }
        if(missing(factor)){
          errorCondition('Only the vector "y" is provided; please also provide "x" and "factor" vectors or provide two lm outputs')
        }
        else{
          errorCondition('Only the vectors "y" and "factor" are provided; please also provide "x" vector or provide two lm outputs')
        }
      }
      else{
        if(missing(factor)){
          errorCondition('Only the vectors "x" and "y" are provided; please also provide "factor" vector or provide two lm outputs')
        }
        if(missing(y)){
          errorCondition('Only the vectors "x" and "factor" are provided; please also provide "y" vector or provide two lm outputs')
        }
        else{
          errorCondition('Only the vector "x" is provided; please also provide "y" and "factor" vectors or provide two lm outputs')
        }
      }
      group_names<-levels(factor(factor))
      object_1<-summary(lm(y[factor==group_names[1]]~x[factor==group_names[1]]))
      object_2<-summary(lm(y[factor==group_names[2]]~x[factor==group_names[2]]))
    }
    if((missing(x)|missing(y)|missing(x2)|missing(y2))==FALSE){
      object_1<-summary(lm(y~x))
      object_2<-summary(lm(y2~x2))
    }
  }
  else{
    if((missing(object_1)==TRUE&missing(object_2)==FALSE)|(missing(object_1)==FALSE)&missing(object_2)==TRUE){
      errorCondition("A single lm output is provided, please provide two of them or provide the original data (x and y data together with a factorial vector distinguishing two groups)")
    }
    if(missing(object_1)==FALSE&missing(object_2)==FALSE){
      if(attributes(object_1)$class=="lm"){
        object_1<-summary(object_1)
      }
      if(attributes(object_2)$class=="lm"){
        object_2<-summary(object_2)
      }
    }
  }
  if(parameter=="slope"){
    parameter<-2
  }
  if(parameter=="intercept"){
    parameter<-1
  }
  parameter_1<-object_1$coefficients[parameter,1]
  parameter_2<-object_2$coefficients[parameter,1]
  se_parameter_1<-object_1$coefficients[parameter,2]
  se_parameter_2<-object_2$coefficients[parameter,2]
  n_1<-length(object_1$residuals)
  n_2<-length(object_2$residuals)
  n<-n_1+n_2
  num<-abs(parameter_1-parameter_2)
  denom<-sqrt(se_parameter_1^2+se_parameter_2^2)
  t<-num/denom
  df<-n-4
  p<-2*pt(-t,df)
  results<-list(num,denom,t,df,p)
  names(results)<-c("Estimate","Std. Error","t value","df","p-value")
  return(results)
}
