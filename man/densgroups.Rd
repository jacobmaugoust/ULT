% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/densgroups.R
\name{densgroups}
\alias{densgroups}
\title{Multi-subgroup density lines}
\usage{
densgroups(x,g,cols,g.levels.order,plot.opt,density.opt,curve.opt,legend=TRUE,legend.opt)
}
\arguments{
\item{x}{The initial variable}

\item{g}{The groups variable, must be of same length than \code{x}}

\item{cols}{The colors of the different variable lines. By default using the first colors of the [ULT:contrasting.palette] function}

\item{g.levels.order}{Optional. The order of the levels for the variable \code{g} if different from the default one}

\item{plot.opt}{Options (as a list) to be passed to the initial \code{plot} function, that sets the global plot frame. By default set to \code{list(xlab=NA,ylab=NA)}}

\item{density.opt}{Options (as a list) to be passed to the \code{density} function used to compute the density.}

\item{curve.opt}{Options (as a list) to be passed to the \code{lines} function used to draw the density lines. By default set to \code{list(lwd=2)}}

\item{legend}{Logical indicating whether a legend for each subgroup has to be drawn}

\item{legend.opt}{Options (as a list) to be passed to the \code{legennd} function if a legend is called. By default set to \code{list("x"="topleft",lwd=2,bty="n",legend=levels(g))}}
}
\description{
This function plots density lines of subgroups for a given variable
}
\examples{
set.seed(1)
means<-runif(5,-20,20)
sds<-runif(5,0.5,3)
x<-c()
for(i in 1:5){x<-c(x,rnorm(100,means[i],sds[i]))}
g<-as.factor(rep(letters[1:5],each=100))
densgroups(x,g) # Default plot
densgroups(x,g,legend.opt=list(x="topleft")) # Default plot changing location of legend
densgroups(x,g, # With some customization
           cols=c("red","blue","brown","pink","cyan"),
           plot.opt=list(xlab="value",ylab=NA,yaxt="n"),
           curve.opt=list(lwd=3,lty=2),
           legend.opt=list("x"="topright",lwd=2,lty=2,bty="n",legend=paste0(levels(g),", mu=",round(means,2),", sd=",round(sds,2))))

}
