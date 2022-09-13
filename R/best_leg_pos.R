#' @title Find best legend position in a plot
#' @description Quick function to automatically add a legend in the plot quarter which contains the less points
#' @param x The x data or the xy dataset
#' @param y Optional. The y data
#' @export
#' @keywords internal

best.leg.pos<-function(x,y){
  if(missing(y)&&ncol(x)==2){
    y<-x[,2]
    x<-x[,1]
  }
  xlim<-par("usr")[1:2]
  ylim<-par("usr")[3:4]
  res<-c("topleft","topright","bottomleft","bottomright")[which.min(c(
    length(which(x<mean(xlim)&y>mean(ylim))),
    length(which(x>mean(xlim)&y>mean(ylim))),
    length(which(x<mean(xlim)&y<mean(ylim))),
    length(which(x>mean(xlim)&y<mean(ylim)))))]
  return(res)
}
