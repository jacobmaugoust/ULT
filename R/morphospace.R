#' @title
#' Function to create a "morphospace"
#'
#' @description
#' This functions performs a "morphospace" of a dataset of two vectors: it connects extremal points to plot the largest polygon possible. This is a rather generic function and can be applied on datasets no matter of which are the vectors. The dataset just needs two vectors: firstly, plot the two vectors against each other (second vector in y-axis; first vector in x-axis), then run morphospace to superimpose a polygon to the plot.
#'
#' @details
#' The function firstly identify the most extreme point on the x-axis, which is the starting point of the polygon.
#' Secondly, the function finds the second point of the polygon which is defined as having the steepest slope with the first point among all other points.
#' Thirdly, the function iteratively calculates the angle (using the same code than the function Angle in the package LearnGeom) between the slopes preceeding and succeeding the i-point (starting from the second defined point): the slope preceeding the i-point is unique (it is defined by the actual i point and the preceeding i-1 point), each slope succeeding the i-point is considered (with all other points). The largest angle is retained together with the associated i+1 point.
#' At the end, the function creates another data.frame of two vectors which each point surrounding the whole scatterplot; the function creates the polygon based on these points.
#'
#' @examples
#'
#' vec1_test<-sample.int(15,replace=TRUE)
#' vec2_test<-sample.int(15,replace=TRUE)
#' df_test<-data.frame(vec1_test,vec2_test)
#' plot(df_test)
#' morphospace(df_test)
#'
#' @param data Data frame with two vectors
#' @param plot.function Type of function used to plot the polygon: "points" or "polygon". Default is "polygon".
#' @param ... graphical arguments, depend of the plot.function choosed
#'
#' @importFrom graphics polygon
#' @importFrom stats na.omit
#'
#' @export
morphospace<-function(data,plot.function,...){
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

  if(missing(plot.function)){
    plot.function<-"points"
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
  }
}
