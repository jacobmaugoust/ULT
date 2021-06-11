#' @title
#' Function to create a "morphospace"
#'
#' @description
#' This functions performs a "morphospace" of a dataset of two vectors: it connects extremal points to plot the largest polygon possible. This is a rather generic function and can be applied on datasets no matter of which are the vectors. The dataset just needs two vectors: firstly, plot the two vectors against each other (second vector in y-axis; first vector in x-axis), then run morphospace to superimpose a polygon to the plot.
#'
#' @details
#' The function firstly identify the most extreme point on the x-axis, which is the starting point of the polygon.
#' Secondly, the function finds the second point of the polygon which is defined as having the steepest slope with the first point among all other points.
#' Thirdly, the function iteratively calculates the angle (using the same code than the function \code{Angle} in the package \code{LearnGeom}) between the slopes preceeding and succeeding the i-point (starting from the second defined point): the slope preceeding the i-point is unique (it is defined by the actual i point and the preceeding i-1 point), each slope succeeding the i-point is considered (with all other points). The largest angle is retained together with the associated i+1 point.
#' At the end, the function creates another data.frame of two vectors which each point surrounding the whole scatterplot; the function creates the polygon based on these points.
#'
#' @seealso
#' The smoothing of the morphospace is adapted from the StackExchange page \href{https://stackoverflow.com/questions/57305921/rasterizing-coordinates-for-an-irregular-polygon-changes-original-shape}{Rasterizing coordinates for an irregular polygon changes original shape}.
#'
#' @examples
#' # For a simple morphospace
#' vec1_test<-sample.int(15,replace=TRUE)
#' vec2_test<-sample.int(15,replace=TRUE)
#' plot(vec1_test,vec2_test)
#' morphospace(vec1_test,vec2_test) # For a "raw" morphospace
#' morphospace(vec1_test,vec2_test,smoothing.method="spline") # For an example of smoothed (here, using splines) morphospace
#' # For a morphospace with several subgroups
#' vec1_test<-c(runif(20,0,10),runif(20,10,20),runif(20,20,30))
#' vec2_test<-sample.int(60,replace=TRUE)
#' groups_test<-as.factor(c(rep("a",20),rep("b",20),rep("c",20)))
#' # For a "raw" morphospace
#' plot(vec1_test,vec2_test,type="n")
#' for (i in 1:length(levels(groups_test))){
#'   morphospace(vec1_test[groups_test==levels(groups_test)[i]],vec2_test[groups_test==levels(groups_test)[i]],pch=21,col=c("red","green3","blue")[i],bg=c("red","green3","blue")[i],plot.type="full")
#' }
#' # For an example of smoothed (here, using splines) morphospace
#' # 1st defining the points of the splines to adjust the size of the plot
#' temp<-matrix(nrow=0,ncol=2,NA)
#' for (i in 1:length(levels(groups_test))){
#'   temp<-rbind(temp,morphospace(vec1_test[groups_test==levels(groups_test)[i]],vec2_test[groups_test==levels(groups_test)[i]],smoothing.method = "spline",output=NA))
#' }
#' # 2nd plotting the data and the morphospaces
#' plot(vec1_test,vec2_test,type="n",xlim=c(min(temp[,1]),max(temp[,1])),ylim=c(min(temp[,2]),max(temp[,2])))
#' for (i in 1:length(levels(groups_test))){ # For an example of smoothed (here, using splines) morphospace
#'   morphospace(vec1_test[groups_test==levels(groups_test)[i]],vec2_test[groups_test==levels(groups_test)[i]],pch=21,col=c("red","green3","blue")[i],bg=c("red","green3","blue")[i],plot.type="full",smoothing.method = "spline")
#' }
#'
#' @param x Either the \code{x} values of the dataset or the whole dataset (in \code{matrix}, \code{data.frame} or \code{function} format).
#' @param y Has to be provided if and only if \code{x} is a single vector.
#' @param plot.function Type of function used to plot the polygon: \code{"points"} or \code{"polygon"}. Default is \code{"points"}.
#' @param plot.type Type of function used to plot the original data if no data has been plotted. By default \code{"lines"}, meaning that no points are plotted; otherwise, \code{"points"} plots the data as points (as in the \code{plot} function).
#' @param output Type of output, either the plot of the current dataset (\code{output="plot"}, the default) or the data allowing to plot the morphospace (\code{output=NA} or any other value) in a little easier way than using \code{chull} function.
#' @param smoothing.method The smoothing method to use if the user wants a smoothed morphospace. See \code{help(smoothr::smooth)} for more details.
#' @param smoothing.param The smoothing parameters to use if the user wants a smoothed morphospace. See \code{help(smoothr::smooth)} for more details. Must be a list with named elements, the names being the names of the parameters.
#' @param plot.new If there has to be a new plot or if morphospace adds to a current plot. By default, it adds to the previous plot.
#' @param ... graphical arguments, depend of the \code{plot.function} choosed
#'
#' @importFrom graphics polygon
#' @importFrom stats na.omit
#' @importFrom rlang is_formula
#' @importFrom sf st_cast
#' @importFrom sf st_multipoint
#' @importFrom smoothr smooth
#' @importFrom stats get_all_vars
#'
#' @export
morphospace<-function(x,y,plot.function,plot.type,output,smoothing.method=NA,smoothing.param=NULL,plot.new=FALSE,...){
  if(missing(y)){
    if(is_formula(x)){
      data<-data.frame(get_all_vars(x))
    }
    else{
      data<-as.data.frame(x)
    }
  }
  else{
    data<-data.frame(x,y)
  }
  data<-na.omit(data)
  points<-data.frame(c(NA),c(NA))
  temp<-data[data[,1]==min(data[,1]),]
  if(is.null(dim(temp))==TRUE){
    points[1,]<-temp
  } else{
    points[1,]<-temp[temp[,2]==min(temp[,2]),]
  }
  pente<-c()
  pentetemp<-c()
  for (i in 1:length(data[,1])){
    if(data[i,1]!=points[1,1]|data[i,2]!=points[1,2]){
      if(is.null(pente)){
        pente<-(data[i,2]-points[1,2])/(data[i,1]-points[1,1])
        points[2,]<-data[i,]
      }
      else{
        pentetemp<-(data[i,2]-points[1,2])/(data[i,1]-points[1,1])
        if(pentetemp>pente){
          pente<-pentetemp
          points[2,]<-data[i,]
        }
      }
    }
  }

  stop<-FALSE
  for (j in 2:length(data[,1])){
    if(stop==TRUE){break}
    else{
      vector1<-c()
      vector2<-c()
      num<-c()
      den<-c()
      angle1<-c()
      first<-c()
      angletemp<-c()
      mark<-FALSE
      for (i in 1:length(data[,1])){
        if(stop==TRUE){break}
        else{
          if(stop==FALSE){
            if(i==length(data[,1])&mark==TRUE){
              stop<-TRUE
              points<-points[-(j+1),]
              break
            }
            if((data[i,1]!=points[j,1]|data[i,2]!=points[j,2])&(data[i,1]!=points[j-1,1]|data[i,2]!=points[j-1,2])){
              if(is.null(angle1)==TRUE){
                vector1<-c(points[j-1,1]-points[j,1],points[j-1,2]-points[j,2])
                vector2<-c(data[i,1]-points[j,1],data[i,2]-points[j,2])
                num<-(vector1[1]*vector2[1]+vector1[2]*vector2[2])
                den<-sqrt(vector1[1]^2+vector1[2]^2)*sqrt(vector2[1]^2+vector2[2]^2)
                angle1<-(360*(acos(num/den)))/(2*pi)
                points[j+1,]<-data[i,]
                first<-i
              }
              else{
                vector2<-c(data[i,1]-points[j,1],data[i,2]-points[j,2])
                num<-(vector1[1]*vector2[1]+vector1[2]*vector2[2])
                den<-sqrt(vector1[1]^2+vector1[2]^2)*sqrt(vector2[1]^2+vector2[2]^2)
                angletemp<-(360*(acos(num/den)))/(2*pi)
                if(angletemp>angle1){
                  if(data[i,1]==points[1,1]&data[i,2]==points[1,2]){
                    mark<-TRUE
                    angle1<-angletemp
                    points[j+1,]<-data[i,]
                  }
                  else{
                    if(mark==TRUE){
                      mark<-FALSE
                      angle1<-angletemp
                      points[j+1,]<-data[i,]
                    }
                    else{
                      angle1<-angletemp
                      points[j+1,]<-data[i,]
                    }
                  }
                  if(i==length(data[,1])&mark==TRUE){
                    stop<-TRUE
                    points<-points[-(j+1),]
                    break
                  }
                }
                else{
                  if(i==length(data[,1])&points[j+1,1]==data[as.numeric(first),1]&data[as.numeric(first),1]==points[1,1]){
                    stop<-TRUE
                    points<-points[-(j+1),]
                    break
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if(is.na(smoothing.method)==FALSE){
    if(any(smoothing.method==c("chaikin", "ksmooth", "spline", "densify"))){
      points<-as.data.frame(matrix(ncol=2,data=unlist(do.call(smooth,c(list(st_cast(st_multipoint(as.matrix(points)),"POLYGON"), method = smoothing.method),smoothing.param)),recursive=FALSE),byrow=F))
    }
    else{
      errorCondition("The chosen smoothing method is not one of those permitted by the smoothr::smooth function; please read the help page")
    }
  }

  if(missing(output)){
    output<-"plot"
  }
  if(output!="plot"|is.na(output)){
    return(points)
  }

  if(missing(plot.type)){
    plot.type<-"lines"
  }

  if(missing(plot.function)){
    plot.function<-"points"
  }

  if(plot.new==TRUE|is.na(smoothing.method)==FALSE){
    plot(data,xlim=c(min(points[,1]),max(points[,1])),ylim=c(min(points[,2]),max(points[,2])),type="n")
  }

  if(plot.function=="polygon"){
    polygon(points,...)
  }
  if(plot.function=="points"){
    for (i in 1:length(points[,1])){
      if(i!=length(points[,1])){
        points(c(points[i,1],points[i+1,1]),c(points[i,2],points[i+1,2]),type="l",...)
      }
      else{
        points(c(points[i,1],points[1,1]),c(points[i,2],points[1,2]),type="l",...)
      }
    }
    if(plot.type!="lines"){
      points(data[,1],data[,2],...)
    }
  }
}
