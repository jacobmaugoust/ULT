#' @title Get a layout grid matrix according to desired number of plots
#'
#' @description This function outputs a matrix to be used with the \link[graphics]{layout} function according to a desired number of plots and to the proportionality between width and height of the plot region
#'
#' @param x The desired number of plots.
#' @param prop A numeric value between 0 and 1. A proportionality factor for distributing plot more or less along plot area width (low value) or height (high value). By default to 0.5, a balanced distribution.
#'
#' @examples
#' # If wanting to do six plots
#' layout.matrix<-get.grid(6)
#' layout(layout.matrix)
#' for(i in 1:6){
#'   plot(1:10,1:10,pch=21,col=NA,bg=ULT::contrasting.palette()[i])
#' }
#' # If wanting to do four plots that are all very tall
#' layout.matrix<-get.grid(4,prop=0)
#' layout(layout.matrix)
#' for(i in 1:4){
#'   boxplot(rnorm(100,i,1))
#' }
#' # If wanting to do four plots that are all very wide
#' layout.matrix<-get.grid(4,prop=1)
#' layout(layout.matrix)
#' for(i in 1:4){
#'   plot(density(rnorm(100,i,1)))
#' }
#' @export

get.grid<-function(x,prop=0.5){
  mults<-unlist(sapply(c(1:x),function(y){if(x/y==round(x/y)){y}}))
  nR<-mults[which.min(abs(prop-sapply(mults,function(y){y/x})))][1]
  nC<-x/nR
  matrix(c(1:x),nrow=nR,ncol=nC,byrow=TRUE)
}
