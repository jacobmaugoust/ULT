#' @title Extended Multivariate Kruskal-Wallis test with missing data
#'
#' @description
#' This function computes an "extended" multivariate Kruskal-Wallis test for n numeric variables (which can contain NA's) relative to one factorial variable (that subsets the dataset in groups)
#'
#' @details
#' "Likelihood-based" and "permutation-based" multivariate Kruskal-Wallis tests are computed: in large samples, the distribution of the test statistic approximates that of the khi?, but in smaller samples, a more accurate p-value is obtained by computing an "empirical" distribution of the test statistic by doing a Monte-Carlo sampling with permutations.
#' Firstly, the "multivariate Kruskal-Wallis test with missing data" is computed and are the first half of the outputs; they are the "likelihood-based" test results (see documentation of \code{multkw.m} for more details).
#' Thus, a Monte-Carlo sampling with permutations (by randomly assigning individuals to groups) is computed, and the second half of the outputs are the proportions of results that exceeds the previously observed results (with the likelihood-based test).
#'
#' @seealso
#' See chapter 2.2.2 and 4.2 of the \href{http://d-scholarship.pitt.edu/19411/1/Fanyin_ETD_draft_08-06-2013.pdf}{PhD manuscript of Fanyin He} and 'Methodology' of \insertCite{He.etal.2017;textual}{ULT} for more details.
#'
#' @param group The factorial variable that subsets the dataset in groups. Can be a character vector, a factorial vector or an integer/numeric vector.
#' @param y The dataset of n numeric(or integer) variables.
#' @param r Optional. The missing data pattern to be applied. If dataset has \code{NA} and if the missing data pattern is the distribution of the \code{NA}'s in the dataset, \code{r} is optional and is automatically computed.
#' @param weight Optional. The weighting scheme to be used to compute the final value of the test statistic. As test statistics are calculated for each pattern of missingness, there are as statistics as patterns. The final test statistic can thus be the arithmetic mean of each statistic (\code{weight="equal"}) or the ponderated mean of each statistic relative to the proportion of each missing pattern (\code{weight="prop"}).
#' @param nmc Number of Monte-Carlo permutations to do.
#' @param print Whether the test should be printed (\code{TRUE}, the default) or not (e.g., to be stored in an object)
#'
#' @importFrom stats aggregate pchisq
#' @import lattice
#' @import Matrix
#' @importFrom Rdpack reprompt
#' @importFrom Rdpack insert_all_ref
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
#' multkw.perm(y=datamkw,airquality$Month,nmc=100)
#' multkw.perm(y=datamkw,airquality$Month,nmc=10000)
#'
#' @export

multkw.perm<-function(nmc,group,y,r,weight,print=TRUE){
  group.var.name<-deparse(substitute(group))
  y.var.name<-deparse(substitute(y))
  if(missing(weight)){weight<-"prop"}
  # count missing patterns #
  g.order<-group
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

  W2.m.perm<-rep(0,nmc)
  W2.c.perm<-rep(0,nmc)

  stats0<-invisible(multkw.m(group,y,r,weight,print=FALSE))
  W2.c<-stats0$W2.c
  nu.c<-stats0$nu.c
  p.multkw.c.chi2<-stats0$p.multkw.c.chi2
  W2.m<-stats0$W2.m
  nu.m<-stats0$nu.m
  p.multkw.m.chi2<-stats0$p.multkw.m.chi2
  pattern.number<-stats0$pattern.number

  for (i in 1:nmc){
    i.st<-1
    group.perm<-rep(0,n)
    group.perm<-sample(group,size=n)
    stats<-invisible(multkw.m(group.perm,y,r,weight,print=FALSE))
    W2.c.perm[i]<-stats$W2.c
    W2.m.perm[i]<-stats$W2.m
  }
  p.multkw.m.perm<-sum(W2.m<W2.m.perm)/nmc
  p.multkw.c.perm<-sum(W2.c<W2.c.perm)/nmc
  multkw.perm.results<-list(y=y.var.name,group=group.var.name,weight=weight,nmc=nmc,pattern.number=pattern.number,W2.c=W2.c,nu.c=nu.c,p.multkw.c.chi2=p.multkw.c.chi2,W2.m=W2.m,nu.m=nu.m,p.multkw.m.chi2=p.multkw.m.chi2,p.multkw.c.perm=p.multkw.c.perm,p.multkw.m.perm=p.multkw.m.perm)
  class(multkw.perm.results)<-"multkw.output"
  if(print){print(multkw.perm.results)}
  else{return(multkw.perm.results)}
}
