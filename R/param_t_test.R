#' @title
#' T-test using describing parameters of data
#'
#' @description
#' Performs one and two sample t-tests using describing parameters (i.e., mean, variance, length) of data
#'
#' @seealso
#' \code{\link[stats]{t.test}}
#'
#' @details
#' \code{alternative = "greater"} is the alternative that \code{x} has a larger mean than \code{y}. For the one-sample case: that the mean is positive.
#' If \code{paired} is \code{TRUE}, an error is returned because a paired t-test can only be done using original data, and the user should perform a t-test using the \code{\link[stats]{t.test}} function
#'
#' @return
#' As for \code{\link[stats]{t.test}}, a list with class \code{"htest"} containing the following components:
#' \describe{
#' \item{\code{statistic}}{the value of the t-statistic}
#' \item{\code{parameter}}{the degrees of freedom of the t-statistic}
#' \item{\code{p.value}}{the p-value for the test}
#' \item{\code{conf.int}}{a confidence interval for the mean appropriate to the specified alternative hypothesis}
#' \item{\code{estimate}}{the (inputted) estimated mean or means depending on whether it was a one-sample or a two-sample test}
#' \item{\code{null.value}}{the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample or a two-sample test}
#' \item{\code{stderr}}{the standard error of the mean/difference, used as denominator in the t-statistic formula}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis}
#' \item{\code{method}}{a character string indicating what type of t-test was performed}
#' \item{\code{data.name}}{a character string giving the (inputted) name(s) of the data}
#' }
#'
#' @export param.t.test
#'
#' @examples
#' ## Define data
#' set.seed(1)
#' A<-runif(10,0,100)
#' B<-runif(10,50,150)
#'
#' mean_x<-mean(A)
#' var_x<-var(A)
#' n_x<-length(A)
#' mean_y<-mean(B)
#' var_y<-var(B)
#' n_y<-length(B)
#'
#' ## One sample t-test
#' t.test(A)
#' param.t.test(mean_x=mean_x,var_x=var_x,n_x=n_x)
#'
#' ## Two sample t-test with equal variances
#' t.test(A,B,var.equal=TRUE)
#' param.t.test(mean_x,var_x,n_x,mean_y,var_y,n_y,var.equal=TRUE)
#'
#' ## Two sample t-test with unequal variances
#' t.test(A,B)
#' param.t.test(mean_x,var_x,n_x,mean_y,var_y,n_y)

param.t.test<-function(mean_x,var_x,n_x,mean_y=NA,var_y=NA,n_y=NA,mu=0,paired=FALSE,var.equal=FALSE,alternative=c("two.sided","less","greater"),conf.level=0.95,name_x="x",name_y="y"){
  if(length(unique(c(is.na(mean_y),is.na(var_y),is.na(n_y))))==2){
    stop("You did not specify all three terms needed for a two-sampled t-test.")
  }
  else{
    if(all(is.na(mean_y),is.na(var_y),is.na(n_y))){
      num<-mean_x-mu
      denom<-sqrt(var_x)/sqrt(n_x)
      df<-n_x-1
      method<-"One Sample t-test"
      estimate<-setNames(mean_x,"mean of x")
      data.name<-name_x
    }
    else{
      num<-mean_x-mean_y
      estimate<-c(setNames(mean_x,"mean of x"),setNames(mean_y, "mean of y"))
      data.name<-paste0(name_x," and ",name_y)
      if(paired==TRUE){
        stop("Paired t-test is only possible with original data, use t.test function")
      }
      else{
        if(var.equal==TRUE){
          denom<-sqrt(((n_x-1)*var_x+(n_y-1)*var_y)/(n_x+n_y-2))*sqrt(1/n_x+1/n_y)
          df<-n_x+n_y-2
          method<-"Two Sample t-test"
        }
        else{
          denom<-sqrt(var_x/n_x+var_y/n_y)
          df<-(var_x/n_x+var_y/n_y)^2/(((var_x/n_x)^2/(n_x-1))+((var_y/n_y)^2/(n_y-1)))
          method<-"Welch Two Sample t-test"
        }
      }
    }
    stat<-num/denom
    if(length(alternative)>1){alternative<-"two.sided"}
    if(alternative=="two.sided"){
      p<-2*pt(-abs(stat),df)
      alpha <- 1 - conf.level
      cint <- qt(1 - alpha/2, df)
      cint <- stat + c(-cint, cint)
    }
    else{
      if(alternative=="greater"){
        p<-pt(stat,df,lower.tail = FALSE)
        cint <- c(stat - qt(conf.level, df), Inf)
      }
      else{
        p<-pt(stat,df)
        cint <- c(-Inf, stat + qt(conf.level, df))
      }
    }
    cint <- mu + cint * denom

    names(stat)<-"t"
    names(df)<-"df"
    conf.int<-cint
    attr(cint, "conf.level") <- conf.level
    names(mu)<-"difference in means"

    results<-list(statistic=stat,
                  parameter=df,
                  p.value=p,
                  conf.int=conf.int,
                  estimate=estimate,
                  null.value=mu,
                  stderr=denom,
                  alternative=alternative,
                  method=method,
                  data.name=data.name)
    class(results)<-"htest"
    results
  }
}
