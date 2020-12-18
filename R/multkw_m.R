#' @title Multivariate Kruskal-Wallis test with missing data
#'
#' @description
#' This function computes a multivariate Kruskal-Wallis test for n numeric variables (which can contain NA's) relative to one factorial variable (that subsets the dataset in groups)
#'
#' @details
#' A "likelihood-based" multivariate Kruskal-Wallis test is computed ; in large samples, the test statistic is approximately khi² distributed.
#' A first "classic" multivariate Kruskal-Wallis test is computed on "complete" data (i.e. removing the rows with at least one missing value).
#' A second test is computed and include missing values: the test is computed for each "missing pattern" (i.e. missing pattern 1 = no missing data ; missing pattern 2 = missing data only in the first variable, etc) and a general test statistic is thus obtained from the "partial" test statistics. See also option "weight".
#' Finally, outputs allow to compare results with complete data only and with missing data.
#' As the test statistic is approximately khi² distributed (in large samples), p-values are based on khi² distributions.
#' Degrees of freedom ar not the same for the "complete" data test and for the "missing" data test, see the See Also section.
#'
#' @seealso
#' See chapter 2.2.2 and 4.2 of the \href{http://d-scholarship.pitt.edu/19411/1/Fanyin_ETD_draft_08-06-2013.pdf}{PhD manuscript of Fanyin He} and 'Methodology' of \insertCite{He.etal.2017;textual}{ULT} for more details.
#'
#' @param group The factorial variable that subsets the dataset in groups. Can be a character vector, a factorial vector or an integer/numeric vector.
#' @param y The dataset of n numeric(or integer) variables.
#' @param r Optional. The missing data pattern to be applied. If dataset has \code{NA} and if the missing data pattern is the distribution of the \code{NA}'s in the dataset, \code{r} is optional and is automatically computed.
#' @param weight Optional. The weighting scheme to be used to compute the final value of the test statistic. As test statistics are calculated for each pattern of missingness, there are as statistics as patterns. The final test statistic can thus be the arithmetic mean of each statistic (\code{weight="equal"}) or the ponderated mean of each statistic relative to the proportion of each missing pattern (\code{weight="prop"}).
#'
#' @importFrom stats aggregate pchisq
#' @import lattice
#' @import Matrix
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{He.etal.2017}{ULT}
#'
#' @return
#' Returns a list of results of the various multivariate Kruskal-Wallis tests that have been computed.
#' The results are the test statistics (W2), the degrees of freedom (df) and the p-value of the test statistic.
#' These three results are given for (1) a "classical" multivariate Kruskal-Wallis test, i.e. on data without missing values; each test statistic is thus followed by a .c for "complete" and (2) a global multivariate Kruskal-Wallis test that takes into account missing values (see details); each test statistic is thus followed by a .m for "missing".
#'
#' @examples
#' data(airquality)
#' datamkw<-airquality[,1:4]
#' multkw(y=datamkw,airquality$Month)
#' multkw.m(y=datamkw,airquality$Month)
#'
#' @export

multkw.m<-function(group,y,r,weight){
  group.var.name<-deparse(substitute(group))
  y.var.name<-deparse(substitute(y))
  if(missing(weight)){weight<-"prop"}
  # count missing patterns #
  g.order<-group
  groupls<-unique(group)
  g<-length(groupls) # number of groups
  p<-dim(y)[2]
  y.order<-y
  if(missing(r)){
    dim_y<-dim(y)
    r<-matrix(as.numeric(is.na(y)),nrow=dim_y[1],ncol=dim_y[2])
  }
  r.order<-r
  for (i in 1:p){
    oo<-order(r.order[,i])
    g.order<-g.order[oo]
    y.order<-y.order[oo,]
    r.order<-r.order[oo,]
  }
  J<-nrow(unique(r.order,MARGIN=1)) # number of missing patterns
  D<-data.frame(r.order)
  n<-length(group)
  ones<-rep(1,n)
  mc<-aggregate(ones,by=as.list(D),FUN=sum) # counts of each missing pattern
  mi<-mc$x
  pi<-p-rowSums(mc[,1:p])

  # get W^2_J #
  W2<-rep(0,J)
  W2.c<-0
  i.st<-1
  for (j in 1:J){
    i.end<-i.st+mi[j]-1
    gg<-g.order[i.st:i.end]
    yy<-y.order[i.st:i.end,]
    ii<-mc[j,1:p]==FALSE
    if(sum(as.numeric(ii))>0){
      yy1<-as.matrix(yy[,ii])
      if(mi[j]>pi[j]){W2[j]<-multkw(gg,yy1)$test.statistic} # if mi[j]>p needs to dig more
    }
    if(prod(as.numeric(ii))==1){W2.c<-W2[j]}
    i.st<-i.end+1
  }
  if(weight=="prop"){tj<-mi/sum(mi)}
  else{tj<-1/J}
  W2.m<-sum(tj*W2)
  denom_nu.m_i<-c()
  for (i in 1:length(pi)){
    if(pi[i]>0){denom_nu.m_i[i]<-(tj[i]*W2.m)^2/pi[i]/(g-1)}
  }
  nu.m<-(W2.m)^2/sum(denom_nu.m_i)
  nu.c<-p*(g-1)
  p.multkw.m.chi2<-pchisq(W2.m,nu.m,lower.tail=FALSE)
  p.multkw.c.chi2<-pchisq(W2.c,nu.c,lower.tail=FALSE)
  multkw.m.results<-list(y=y.var.name,group=group.var.name,weight=weight,pattern.number=length(pi),W2.c=W2.c,nu.c=nu.c,p.multkw.c.chi2=p.multkw.c.chi2,W2.m=W2.m,nu.m=nu.m,p.multkw.m.chi2=p.multkw.m.chi2)
  class(multkw.m.results)<-"multkw.output"
  return(multkw.m.results)
}
