library(lattice)
library(Matrix)

# 1. mult-KW function

multkw<- function(group,y,simplify=FALSE){
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

  # return stat and p-value #
  returnlist<-list(statistic=W2,d.f.=p*(g-1),p.value=pchisq(W2,p*(g-1),lower.tail = FALSE))
  if(simplify==TRUE){return(W2)}
  else{return(returnlist)}
}

# 2. multkw with missing values

multkw.m<-function(group,y,r,weight){
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
      if(mi[j]>pi[j]){W2[j]<-multkw(gg,yy1,simplify=TRUE)} # if mi[j]>p needs to dig more
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
  return(list(W2.m=W2.m,nu.m=nu.m,p.multkw.m.chi2=p.multkw.m.chi2,W2.c=W2.c,nu.c=nu.c,p.multkw.c.chi2=p.multkw.c.chi2)) # W2.c = W2 for complete cases only ; nu = degrees of freedom following Welch-Satterthwaite equation where W2 is approximately X² distributed ; W2.m = W2 for missing values = for all cases except those where there are only NA's
}

# 3. Monte-Carlo permutations

multkw.perm<-function(nmc,group,y,r,weight){
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
  stats0<-multkw.m(group,y,r,weight)
  W2.m<-stats0$W2.m
  W2.c<-stats0$W2.c
  nu<-stats0$nu
  for (i in 1:nmc){
    i.st<-1
    group.perm<-rep(0,n)
    group.perm<-sample(group,size=n)
    stats<-multkw.m(group.perm,y,r,weight)
    W2.m.perm[i]<-stats$W2.m
    W2.c.perm[i]<-stats$W2.c
  }
  p.multkw.m.perm<-sum(W2.m<W2.m.perm)/nmc
  p.multkw.m.chi2<-pchisq(W2.m,nu,lower.tail=FALSE)
  p.multkw.c.perm<-sum(W2.c<W2.c.perm)/nmc
  p.multkw.c.chi2<-pchisq(W2.c,p*(g-1),lower.tail=FALSE)
  return(list(W2.m=W2.m,p.multkw.m.perm=p.multkw.m.perm,p.multkw.m.chi2=p.multkw.m.chi2,p.multkw.c.perm=p.multkw.c.perm,p.multkw.c.chi2=p.multkw.c.chi2))
}

# 4. Data generation

data.gen<-function(p,g=2,ni,delta,mpcnt){
  n<-ni*g
  X1<-rnorm(ni)
  Y11<-sapply(X1+1,rnorm,n=1,sd=sqrt(2))
  Y12<-sapply(X1,rnorm,n=1,sd=1)
  X2<-rnorm(ni)
  Y21<-sapply(X2+1,rnorm,n=1,sd=sqrt(2))
  Y22<-sapply(X2+delta,rnorm,n=1,sd=1)
  Y1<-matrix(c(Y11,Y12),nrow=ni,ncol=p)
  Y2<-matrix(c(Y21,Y22),nrow=ni,ncol=p)
  y<-rbind(Y1,Y2)
  group<-rep(1:g,each=ni)
  m1<-matrix(rep(c(0,0),times=n*mpcnt[1]),ncol=p,byrow=TRUE)
  m2<-matrix(rep(c(0,1),times=n*mpcnt[2]),ncol=p,byrow=TRUE)
  m3<-matrix(rep(c(1,0),times=n*mpcnt[3]),ncol=p,byrow=TRUE)
  m<-rbind(m1,m2,m3)
  perm<-sample(n)
  r<-m[perm,]
  return(list(group=group,y=y,r=r))
}

data.gen2<-function(p,g=2,ni,delta,mpcnt){
  n<-ni*g
  X1<-rbinom(ni,5,0.5)
  W1<-rbinom(ni,2,0.5)
  Y11<-sapply(X1+1,rpois,n=1)
  Y12<-sapply(X1+2,rpois,n=1)
  X2<-rbinom(ni,5,0.5)
  Y21<-sapply(X2+1,rpois,n=1)
  Y22<-sapply(X2+2+delta,rpois,n=1)
  Y1<-matrix(c(Y11,Y12),nrow=ni,ncol=p)
  Y2<-matrix(c(Y21,Y22),nrow=ni,ncol=p)
  y<-rbind(Y1,Y2)
  group<-rep(1:g,each=ni)
  m1<-matrix(rep(c(0,0),times=n*mpcnt[1]),ncol=p,byrow=TRUE)
  m2<-matrix(rep(c(0,1),times=n*mpcnt[2]),ncol=p,byrow=TRUE)
  m3<-matrix(rep(c(1,0),times=n*mpcnt[3]),ncol=p,byrow=TRUE)
  m<-rbind(m1,m2,m3)
  perm<-sample(n)
  r<-m[perm,]
  return(list(group=group,y=y,r=r))
}

# 5. Simulation

p<-2
g<-2
ni<-50
delta<-2.5
mpcnt<-c(0.2,0.4,0.4)

nsim<-10
psim<-matrix(0,nrow=nsim,ncol=8)
wilks<-rep(0,nsim)
tau<-rep(0,nsim)

nmc<-10
r.comp<-matrix(0,nrow=ni*g,ncol=p)

chi2<-matrix(0,nrow=10000,ncol=3)
chi2[,1]<-rchisq(10000,df=2)
chi2[,2]<-rchisq(10000,df=1)
chi2[,3]<-rchisq(10000,df=1)
ptable1<-rowMeans(chi2)
ptable2<-chi2%*%as.matrix(mpcnt)

Sys.time()
for (i in 1:nsim){
  data1<-data.gen2(p,g,ni,delta,mpcnt)
  group<-data1$group
  y<-data1$y
  r<-data1$r
  fit<-manova(y~group)
  wilks[i]<-(summary(fit,test="Wilks"))$stats[1,2] # value of Wilks' lambda
  ps1<-multkw.perm(nmc,group,y,r,weight="plain")
  ps2<-multkw.perm(nmc,group,y,r,weight="prop")
  ps.comp<-multkw.perm(nmc,group,y,r.comp,weight="prop")
  psim[i,1]<-ps1$p.multkw.m.perm
  psim[i,2]<-mean(ps1$W2.m<ptable1)
  psim[i,3]<-ps1$p.multkw.c.perm
  psim[i,4]<-ps1$p.multkw.c.chi2
  psim[i,5]<-ps2$p.multkw.m.perm
  psim[i,6]<-mean(ps2$W2.m<ptable2)
  psim[i,7]<-ps.comp$p.multkw.m.perm
  psim[i,8]<-ps.comp$p.multkw.m.chi2
}
Sys.time()

pw<-colMeans(psim<0.05)
f2<-mean(1/wilks-1)

fname<-paste("power est n",ni,"delta",delta,"-poi-h",".txt")
write.table(c(f2,pw),file=fname,row.names=TRUE,append=TRUE)
