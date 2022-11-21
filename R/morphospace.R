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
#' morphospace(vec1_test,vec2_test,smoothing.method="spline",plot.points = FALSE) # For an example of smoothed (here, using splines) morphospace
#' morphospace(vec1_test,vec2_test,smoothing.method="spline",plot.points = TRUE) # Same, with the data
#' # For a morphospace with several subgroups
#' vec1_test<-c(runif(20,0,10),runif(20,10,20),runif(20,20,30))
#' vec2_test<-sample.int(60,replace=TRUE)
#' groups_test<-as.factor(c(rep("a",20),rep("b",20),rep("c",20)))
#' # For a "raw" morphospace
#' plot(vec1_test,vec2_test,type="n")
#' morphospace(vec1_test,vec2_test,groups_test,pch=21,col=c("red","green3","blue"),bg=c("red","green3","blue"))
#' # For an example of smoothed (here, using splines) morphospace
#' morphospace(vec1_test,vec2_test,groups_test,pch=21,col=c("red","green3","blue"),bg=c("red","green3","blue"),smoothing.method="spline")
#' # The same morphospace using the polygon function and its options to draw dotted semi-transparent lines
#' morphospace(vec1_test,vec2_test,groups_test,plot.function="polygon",pch=21,col=scales::alpha(c("red","green3","blue"),0.25),density=15,border=c("red","green3","blue"),col.points=c("red","green3","blue"),bg=c("red","green3","blue"),lty=2,plot.new=TRUE)
#' # The same morphospace using the polygon function and its options to draw semi-transparent convex hulls
#' morphospace(vec1_test,vec2_test,groups_test,plot.function="polygon",pch=21,col=scales::alpha(c("red","green3","blue"),0.1),border=c("red","green3","blue"),col.points=c("red","green3","blue"),bg=c("red","green3","blue"),lty=2,plot.new=TRUE)
#'
#' @param x Either the \code{x} values of the dataset or the whole dataset (in \code{matrix}, \code{data.frame} or \code{function} format).
#' @param y Has to be provided if and only if \code{x} is a single vector.
#' @param groups Optional. A factor vector defining subsets (or 'groups') in the dataset to plot separate morphospaces at once.
#' @param plot.function Type of function used to plot the polygon: \code{"points"} or \code{"polygon"}. Default is \code{"points"}.
#' @param plot.points Logical. If the data points have to be plotted. Default is \code{TRUE}.
#' @param output Type of output, either the plot of the current dataset (\code{output="plot"}, the default) or the data allowing to plot the morphospace (\code{output=NA} or any other value) in a little easier way than using \code{chull} function.
#' @param smoothing.method The smoothing method to use if the user wants a smoothed morphospace. See \link[smoothr]{smooth} for more details.
#' @param smoothing.param The smoothing parameters to use if the user wants a smoothed morphospace. See \link[smoothr]{smooth} for more details. Must be a list with named elements, the names being the names of the parameters.
#' @param plot.new If there has to be a new plot or if morphospace adds to a current plot. By default, it adds to the previous plot.
#' @param plot.new.opt The options to be used if a new plot is added (if requested or if drawing smoothed morphospaces). Only works for the general frame of the plot (axis, labels etc). Must be a list of arguments.
#' @param ... graphical arguments, depend of the \code{plot.function} choosed. If \code{plot.function="polygon"} and \code{plot.points=TRUE} and the user wants to specify the \code{col} arguments for both the polygon (i.e., filling color) and the points (i.e., border points color), it is recommended to specify the polygon \code{col} as \code{col} and the points \code{col} as \code{col.points}.
#'
#' @importFrom graphics polygon
#' @importFrom stats na.omit
#' @importFrom rlang is_formula
#' @importFrom sf st_cast
#' @importFrom sf st_multipoint
#' @importFrom smoothr smooth
#' @importFrom stats get_all_vars
#' @importFrom foreach foreach
#'
#' @export
morphospace<-function(x,y,groups,plot.function,plot.points=TRUE,output="plot",smoothing.method=NA,smoothing.param=NULL,plot.new=FALSE,plot.new.opt=NULL,...){
  if(missing(y)){
    if(is_formula(x)){
      orig_data<-data.frame(get_all_vars(x))
    }
    else{
      orig_data<-as.data.frame(x)
    }
  }
  else{
    orig_data<-data.frame(x,y)
  }

  na_data<-FALSE
  if(any(is.na(orig_data))){
    na_data<-TRUE
  }

  if(!missing(groups)){
    groups<-as.factor(groups)
    if(length(groups)!=length(orig_data[,1])&!na_data){
      warning("group vector is not the same length than the quantitative data; groups discarded, please provide a vector of same length")
      groups<-rep(as.factor("one"),length(orig_data[,1]))
    }
  }
  else{
    groups<-rep(as.factor("one"),length(orig_data[,1]))
  }

  if(na_data){
    orig_data<-na.omit(orig_data)
    if(length(groups)!=length(orig_data[,1])){
      groups<-groups[-attr(orig_data,"na.action")]
    }
  }

  repeats<-c()
  for (i in 2:length(orig_data[,1])){
    if(all(orig_data[i,]==orig_data[i-1,])){
      repeats<-c(repeats,i)
    }
  }
  if(length(repeats)>0){
    orig_data<-orig_data[-repeats,]
    groups<-groups[-repeats]
  }

  old_groups<-groups
  oldngroups<-nlevels(old_groups)

  singleton<-c()
  if(any(table(groups)==1)){
    singleton<-names(table(groups))[which(table(groups)==1)]
    warning(paste0("The following groups have not been represented because they consist of a single point: ",ifelse(length(singleton)>1,paste0(singleton,collapse=", "),singleton)))
  }
  if(length(singleton)>0){
    groups_data<-orig_data[-which(groups%in%singleton),]
    groups<-groups[-which(groups%in%singleton)]
    groups<-factor(groups)
  }
  else{
    groups_data<-orig_data
  }

  ngroups<-nlevels(groups)
  groupslevels<-levels(groups)
  points<-list()
  min_x<-c()
  max_x<-c()
  min_y<-c()
  max_y<-c()

  for (i in 1:ngroups){
    points[[i]]<-data.frame(c(NA),c(NA))
    data<-groups_data[groups==groupslevels[i],]
    temp<-data[data[,1]==min(data[,1]),]
    if(is.null(dim(temp))==TRUE){
      points[[i]][1,]<-temp
    } else{
      points[[i]][1,]<-temp[temp[,2]==min(temp[,2]),]
    }
    pente<-c()
    pentetemp<-c()
    for (j in 1:length(data[,1])){
      if(data[j,1]!=points[[i]][1,1]|data[j,2]!=points[[i]][1,2]){
        if(is.null(pente)){
          pente<-(data[j,2]-points[[i]][1,2])/(data[j,1]-points[[i]][1,1])
          points[[i]][2,]<-data[j,]
        }
        else{
          pentetemp<-(data[j,2]-points[[i]][1,2])/(data[j,1]-points[[i]][1,1])
          if(pentetemp>pente){
            pente<-pentetemp
            points[[i]][2,]<-data[j,]
          }
        }
      }
    }

    if(length(data[,1])>2){
      stop<-FALSE
      for (j in 2:length(data[,1])){
        vector1<-c()
        vector2<-c()
        num<-c()
        den<-c()
        angle1<-c()
        first<-c()
        angletemp<-c()
        if(stop){break}else{
          for (k in 1:length(data[,1])){
            if(!all(data[k,]==points[[i]][j,])&!all(data[k,]==points[[i]][j-1,])){
              if(is.null(angle1)){
                vector1<-c(points[[i]][j-1,1]-points[[i]][j,1],points[[i]][j-1,2]-points[[i]][j,2])
                vector2<-c(data[k,1]-points[[i]][j,1],data[k,2]-points[[i]][j,2])
                num<-(vector1[1]*vector2[1]+vector1[2]*vector2[2])
                den<-sqrt(vector1[1]^2+vector1[2]^2)*sqrt(vector2[1]^2+vector2[2]^2)
                angle1<-(360*(acos(num/den)))/(2*pi)
                points[[i]][j+1,]<-data[k,]
              }
              else{
                vector2<-c(data[k,1]-points[[i]][j,1],data[k,2]-points[[i]][j,2])
                num<-(vector1[1]*vector2[1]+vector1[2]*vector2[2])
                den<-sqrt(vector1[1]^2+vector1[2]^2)*sqrt(vector2[1]^2+vector2[2]^2)
                angletemp<-(360*(acos(num/den)))/(2*pi)
                if(angletemp>angle1){
                  angle1<-angletemp
                  points[[i]][j+1,]<-data[k,]
                }
              }
            }
            if(k==length(data[,1])&all(points[[i]][j+1,]==points[[i]][1,])){
              points[[i]]<-points[[i]][-(j+1),]
              stop<-TRUE
              break
            }
          }
        }
      }
    }

    if(is.na(smoothing.method)==FALSE){
      if(any(smoothing.method==c("chaikin", "ksmooth", "spline", "densify"))){
        points[[i]]<-as.data.frame(matrix(ncol=2,data=unlist(do.call(smoothr::smooth,c(list(st_cast(st_multipoint(as.matrix(points[[i]])),"POLYGON"), method = smoothing.method),smoothing.param)),recursive=FALSE),byrow=F))
      }
      else{
        errorCondition("The chosen smoothing method is not one of those permitted by the smoothr::smooth function; please read the help page")
      }
    }

    min_x<-c(min_x,min(points[[i]][,1]))
    max_x<-c(max_x,max(points[[i]][,1]))
    min_y<-c(min_y,min(points[[i]][,2]))
    max_y<-c(max_y,max(points[[i]][,2]))
  }

  names(points)<-as.character(c(1:oldngroups)[levels(old_groups)%in%groupslevels])

  if(output!="plot"|is.na(output)){
    return(points)
  }

  else{
    if(missing(plot.function)){
      plot.function<-"points"
    }

    if(plot.new==TRUE|is.na(smoothing.method)==FALSE){
      if(!"xlim"%in%names(plot.new.opt)){
        plot.new.opt<-c(plot.new.opt,list(xlim=c(min(min_x),max(max_x))))
      }
      if(!"ylim"%in%names(plot.new.opt)){
        plot.new.opt<-c(plot.new.opt,list(ylim=c(min(min_y),max(max_y))))
      }
      if(plot.points==FALSE){
        plot.new.opt<-c(plot.new.opt,list("type"="n"))
      }
      do.call(plot,c(list(orig_data),plot.new.opt))
    }

    args<-list(...)
    if(length(args)>0){
      singles<-c()
      multiples<-c()
      for (i in 1:length(args)){
        if(length(args[[i]])>1){
          multiples<-c(multiples,i)
        }
        else{
          singles<-c(singles,i)
        }
      }
    }

    for (i in 1:oldngroups){
      if(length(args)>0){
        i_args<-c(setNames(foreach(j=multiples)%do%args[[j]][i],names(args)[multiples]),args[singles])
      }
      else{
        i_args<-NULL
      }

      if(levels(old_groups)[i]%in%groupslevels){

      }

      lines_i_args<-c(list("x"=points[[as.character(i)]][,1]),list("y"=points[[as.character(i)]][,2]))
      if(plot.function=="polygon"){
        lines_i_args<-c(lines_i_args,i_args[names(i_args)%in%c("density","angle","border","col","lty","xpd","lend","ljoin","lmitre","fillOddEven")])
        if(plot.points){
          points_i_args<-i_args[!names(i_args)%in%c("density","angle","border","col","lty","xpd","lend","ljoin","lmitre","fillOddEven")]
        }
      }
      else{
        lines_i_args<-c(lines_i_args,i_args,list("type"="l"))
        if(!is.null(lines_i_args$x)){
          lines_i_args$x<-c(lines_i_args$x,points[[as.character(i)]][1,1])
          lines_i_args$y<-c(lines_i_args$y,points[[as.character(i)]][1,2])
        }
        if(plot.points){
          points_i_args<-i_args
        }
      }

      do.call(plot.function,lines_i_args)

      if(plot.points){
        toplot_data<-orig_data[old_groups==levels(old_groups)[i],]
        if(any(names(points_i_args)=="col.points")){
          names(points_i_args)[names(points_i_args)=="col.points"]<-"col"
        }
        do.call("points",c(list("x"=toplot_data[,1]),list("y"=toplot_data[,2]),points_i_args))
      }
    }
  }
}
