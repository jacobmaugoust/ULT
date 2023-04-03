#' @title Swap dimensions of a list
#'
#' @description This function inverts the dimension of a given list; for instance, columns of each elements can become the new listing. All elements of the list must be matrices of same dimensions.
#'
#' @param list The list whose dimensions have to be swapped
#' @param new.dim A numeric vector specifying the re-ordering. Hierarchy of elements is: rows, columns, elements (for i rows, j columns, and k elements, one has \code{list[[k]][i,j]})
#' @param progress.bar A logical indicating whether a progression bar should be printed or not
#'
#' @importFrom stats setNames
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @examples
#' # Create a simple list
#' list<-list("a"=matrix(1:15,ncol=3,nrow=5),"b"=matrix((1:15)^2,ncol=3,nrow=5),"c"=matrix((1:15)/2,ncol=3,nrow=5),"d"=matrix(log(1:15),ncol=3,nrow=5))
#' # Elements are named; set x1, x2.. and y1, y2... as row and column names for each matrix (so that we can ensure the swapping is done properly)
#' for(i in 1:length(list)){dimnames(list[[i]])<-list(paste("x",1:nrow(list[[i]]),sep=""),paste("y",1:ncol(list[[i]]),sep=""))}
#' list
#' # If wanting a list of 3 elements, each being a matrix of five columns and four rows (so that y are elements, x are columns, and a-b-c-d are rows)
#' swap.dim.list(list,new.dim=c(3,1,2)) # The new lines (first element of new.dim) are the previous elements, the new columns are the previous rows, and the new elements are the previous columns
#' # Same if wanting x to be elements, a-b-c-d to be columns, and y to be rows (so former columns become rows, former elements become columns, and former rows become elements)
#' swap.dim.list(list,new.dim=c(2,3,1))
#' @export

swap.dim.list<-function(list,new.dim=c(1,2,3),progress.bar=TRUE){
  lengths<-c(nrow(list[[1]]),ncol(list[[1]]),length(list))
  names<-list(rownames(list[[1]]),colnames(list[[1]]),names(list))
  nR<-lengths[new.dim[1]]
  nC<-lengths[new.dim[2]]
  nL<-lengths[new.dim[3]]
  res<-list(matrix(NA,ncol=nC,nrow=nR,dimnames=list(names[[new.dim[1]]],names[[new.dim[2]]])))
  res<-setNames(rep(res,nL),names[[new.dim[3]]])
  if(progress.bar){
    pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                         max = nR*nC*nL, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 100,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    p<-0
  }
  for(i in 1:nrow(list[[1]])){
    for(j in 1:ncol(list[[1]])){
      for(k in 1:length(list)){
        res[[eval(parse(text=c("i","j","k")[new.dim[3]]))]][eval(parse(text=c("i","j","k")[new.dim[1]])),eval(parse(text=c("i","j","k")[new.dim[2]]))]<-list[[k]][i,j]
        if(progress.bar){
          p<-p+1
          setTxtProgressBar(pb, p)
        }
      }
    }
  }
  if(progress.bar){
    close(pb)
  }
  res
}
