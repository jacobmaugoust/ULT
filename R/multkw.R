#' @title Multivariate Kruskal-Wallis test
#'
#' @description This function computes a multivariate Kruskal-Wallis test for n numeric variables relative to one factorial variable (that subsets the dataset in groups)
#'
#' @details
#' A "standard" multivariate Kruskal-Wallis test is computed, deleting all missing data.
#'
#' @seealso
#' See chapter 2.2.2 and 4.2 of the \href{http://d-scholarship.pitt.edu/19411/1/Fanyin_ETD_draft_08-06-2013.pdf}{PhD manuscript of Fanyin He} and 'Methodology' of \insertCite{He.etal.2017;textual}{ULT} for more details.
#'
#' @param group The factorial variable that subsets the dataset in groups. Can be a character vector, a factorial vector or an integer/numeric vector.
#' @param y The dataset of n numeric(or integer) variables.
#' @param print Whether the test should be printed (\code{TRUE}, the default) or not (e.g., to be stored in an object)
#'
#' @return Output is either a list (with \code{"simplify=FALSE"}) or a vector (with \code{"simplify=TRUE"}) containing the results of the multivariate Kruskal-Wallis test.
#'
#' @importFrom stats aggregate pchisq
#' @import lattice
#' @import Matrix
#'
#' @references
#' \insertRef{He.etal.2017}{ULT}
#'
#' @examples
#' data(airquality)
#' datamkw<-airquality[,1:4]
#' multkw(y=datamkw,airquality$Month)
#'
#' @export

multkw<- function(group,y,print=TRUE){
  group.var.name<-deparse(substitute(group))
  y.var.name<-deparse(substitute(y))
  # sort and rank data by group #
  o<-order(group)
  group<-group[o]
  if(ncol(as.matrix(y))==1){y<-as.matrix(y[o])}
  else{y<-as.matrix(y[o,])}
  n<-length(group)
  p<-dim(y)[2]
  if(dim(y)[1]!=n){return("number of observations not equal to length of group")}
  groupls<-unique(group)
  g<-length(groupls) # number of groups
  groupind<-sapply(groupls,"==",group) # group indicator
  ni<-colSums(groupind) # number of individuals of each group
  r<-apply(y,2,rank) # corresponding rank variable

  # calculation of statistic #
  r.ik<-t(groupind)%*%r*(1/ni) # gxp, mean rank of k-th variate in i-th group
  m<-(n+1)/2 # expected value of rik
  u.ik<-t(r.ik-m)
  U<-as.vector(u.ik)
  V<-1/(n-1)*t(r-m)%*%(r-m) # pooled within-group cov matrix
  Vstar<-bdiag(lapply(1/ni,"*",V))
  W2<-as.numeric(t(U)%*%solve(Vstar)%*%U)
  df<-p*(g-1)
  pv<-pchisq(W2,p*(g-1),lower.tail = FALSE)

  multkw.results<-list(y=y.var.name,group=group.var.name,test.statistic=W2,df=df,p.value=pv)
  class(multkw.results)<-"multkw.output"
  if(print){print(multkw.results)}
  else{return(multkw.results)}
}
