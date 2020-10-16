#' @title T-test of two regression slopes or intercepts
#'
#' @description
#' This function performs a t-test to test whether two regression slopes or two regression intercepts are equal or not.
#' It can also be used to test whether a regression slope or intercept is equal to a given value.
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
#' ## To test a pre-established vs estimated parameter
#' x<-c(sort(runif(100,0,100)))
#' y<-c()
#' slope<-runif(1,-100,100)
#' intercept<-runif(1,-100,100)
#' for (i in 1:length(x)){
#'   y[i]<-x[i]*slope+intercept+runif(1,-abs(slope)/100,abs(slope)/100)
#' }
#' parameter<-"slope"
#' reg.param.t.test(x,y,slope,parameter)
#' reg.param.t.test(y~x,slope,"slope")
#' reg.param.t.test(lm(y~x),slope,parameter) ## This tests if the function retrieves the original defined slope as equal to the slope of
#' reg.param.t.test(lm(y~x),intercept,"intercept") ## This tests if the function retrieves the original defined intercept as equal to the intercept of lm()
#' reg.param.t.test(lm(y~x),0,parameter)[c(1,2,3,5)]==summary(lm(y~x))$coefficients[2,]
#' reg.param.t.test(lm(y~x),0,"intercept")[c(1,2,3,5)]==summary(lm(y~x))$coefficients[1,]
#'
#' ## To compare two pre-established slopes
#' x1<-sort(runif(50,0,20))
#' x2<-sort(runif(50,0,20))
#' slope1<-runif(1,-100,100)
#' intercept<-runif(1,-100,100)
#' slope2<-runif(1,-100,100)
#' for (i in 1:length(x1)){
#'   y1[i]<-x1[i]*slope1+intercept+runif(1,-abs(slope1)/100,abs(slope1)/100)
#'   y2[i]<-x2[i]*slope2+intercept+runif(1,-abs(slope2)/100,abs(slope2)/100)
#' }
#' x12<-c(x1,x2)
#' y12<-c(y1,y2)
#' fac<-c(rep("1",50),rep("2",50))
#' object_1<-summary(lm(y12[fac==levels(factor(fac))[1]]~x12[fac==levels(factor(fac))[1]]))
#' object_2<-summary(lm(y12[fac==levels(factor(fac))[2]]~x12[fac==levels(factor(fac))[2]]))
#' plot(x1,y1,type="l",col="blue",ylim=c(min(min(y1),min(y2)),max(max(y1),max(y2))))
#' points(x2,y2,type="l",col="red")
#' reg.param.t.test(x12,y12,fac,parameter)
#' reg.param.t.test(y12~x12,fac,parameter)
#' reg.param.t.test(y1~x1,y2~x2,parameter)
#' reg.param.t.test(x1,y1,x2,y2,parameter)
#' reg.param.t.test(object_1,object_2,parameter)
#' summary(lm(y12~x12*fac))$coefficients[4,] ## Results do not perfectly fit but are roughly the same, don't know why now...
#'
#' @importFrom rlang is_formula
#' @export

reg.param.t.test<-function(x,y,factor,x2,y2,object_1=NULL,object_2=NULL,delta=NULL,parameter="slope"){

  if(missing(y2)==FALSE){
    if(y2=="slope"|y2=="intercept"){
      parameter<-y2
      y2<-NULL
    }
  }
  if(missing(y2)){y2<-NULL}

  if(missing(x2)==FALSE){
    if((is.numeric(x2)&length(x2)>1)&is.null(y2)){
      y2<-x2
      x2<-NULL
    }
    else if(x2=="slope"|x2=="intercept"){
      parameter<-x2
      x2<-NULL
    }
  }
  if(missing(x2)){x2<-NULL}

  if(missing(factor)==FALSE){
    if(is.numeric(factor)){
      if(length(factor)==1){
        delta<-factor
        factor<-NULL
      }
      else{
        x2<-factor
        factor<-NULL
      }
    }
    else if(length(factor)==1){
      parameter<-factor
      factor<-NULL
    }
  }
  if(missing(factor)){factor<-NULL}

  if(missing(y)==FALSE){
    if(is.atomic(y)){
      if(is.numeric(y)&length(y)==1){
        delta<-y
        y<-NULL
      }
      else if((is.factor(y)|is.character(y))&length(y)>1){
        factor<-y
        y<-NULL
      }
    }
    else if(rlang::is_formula(y)){
      object_2<-summary(lm(y))
      y<-NULL
    }
    else if(is.list(y)&(exists(attributes(y)$class)&(attributes(y)$class=="lm"|attributes(y)$class=="summary.lm"))){
      object_2<-y
      y<-NULL
    }
  }

  if(is.atomic(x)==FALSE){
    if(rlang::is_formula(x)){
      if(is.null(object_2)){
        x_lm<-all.vars(x)
        x<-get(x_lm[2], envir = parent.frame())
        y<-get(x_lm[1], envir = parent.frame())
      }
      else{
        object_1<-summary(lm(x))
        x<-NULL
      }
    }
    else if(is.list(x)&(exists(attributes(x)$class)&(attributes(x)$class=="lm"|attributes(x)$class=="summary.lm"))){
      object_1<-x
      x<-NULL
    }
  }

  dfsup<-0
  if(is.null(object_1)==FALSE&is.null(object_2)&is.null(delta)==FALSE){
    object_2<-list(coefficients=data.frame(c(delta,delta),c(0,0)),residuals=c(0))
    attributes(object_2)$class<-"summary.lm"
    dfsup<--1
  }
  if((is.null(object_1)==FALSE|is.null(object_2)==FALSE)&(is.null(x)==FALSE|is.null(y)==FALSE|is.null(factor)==FALSE|is.null(x2)==FALSE|is.null(y2)==FALSE)==FALSE){
    errorCondition("You provided a lm output and a vector; please choose the method to use and provide either two lm outputs or x and y data together with a factorial vector distinguishing two groups")
  }
  if(is.null(object_1)&is.null(object_2)){
    if(is.null(x2)&is.null(y2)){
      if(is.null(x)){
        if(is.null(y)){
          errorCondition('Only the vector "factor" is provided; please also provide "x" and "y" vectors or provide two lm outputs')
        }
        if(is.null(factor)){
          errorCondition('Only the vector "y" is provided; please also provide "x" and "factor" vectors or provide two lm outputs')
        }
        else{
          errorCondition('Only the vectors "y" and "factor" are provided; please also provide "x" vector or provide two lm outputs')
        }
      }
      else{
        if(is.null(y)==FALSE&(is.null(factor)&is.null(delta))){
          errorCondition('Only the vectors "x" and "y" are provided; please also provide either "factor" vector, a "delta" level, or two lm outputs')
        }
        if(is.null(y)==TRUE&(is.null(factor)==FALSE|is.null(delta)==FALSE)){
          errorCondition('Only the vectors "x" and "factor" or "delta" are provided; please also provide "y" vector or provide two lm outputs')
        }
        if(is.null(y)&(is.null(factor)&is.null(delta))){
          errorCondition('Only the vector "x" is provided; please also provide "y" and "factor" vectors or delta level, or provide two lm outputs')
        }
      }
      if(is.null(factor)==FALSE&is.null(delta)){
        group_names<-levels(factor(factor))
        object_1<-summary(lm(y[factor==group_names[1]]~x[factor==group_names[1]]))
        object_2<-summary(lm(y[factor==group_names[2]]~x[factor==group_names[2]]))
      }
      else{
        object_1<-summary(lm(y~x))
        object_2<-list(coefficients=data.frame(c(delta,delta),c(0,0)),residuals=c(0))
        attributes(object_2)$class<-"summary.lm"
        dfsup<--1
      }
    }
    else if((is.null(x)&is.null(y)&is.null(x2)&is.null(y2))==FALSE){
      object_1<-summary(lm(y~x))
      object_2<-summary(lm(y2~x2))
    }
  }
  else{
    if(is.null(object_1)==FALSE&is.null(delta)==TRUE){
      errorCondition("A single lm output is provided, please provide a comparison value (delta), or provide either two lm output or the original data (x and y data together with a factorial vector distinguishing two groups)")
    }
    if(is.null(object_1)==FALSE){
      if(attributes(object_1)$class=="lm"){
        object_1<-summary(object_1)
      }
    }
    if(is.null(object_1)==FALSE){
      if(attributes(object_2)$class=="lm"){
        object_2<-summary(object_2)
      }
    }
  }
  if(parameter=="slope"){
    parameter<-2
  }
  else if(parameter=="intercept"){
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
  df<-n-4-dfsup
  p<-2*pt(-t,df)
  results<-list(num,denom,t,df,p)
  names(results)<-c("Estimate","Std. Error","t value","df","p-value")
  return(results)
}
