#' @title
#' Get the distance matrix of string data
#'
#' @description
#' This function computes a distance matrix (as does the function \link[stats]{dist}) of character values of a dataset (while \link[stats]{dist} does for numerical data).
#'
#' @details
#' The function counts the number of differences between the rows of a dataset.
#' Please note that this function also works with non-character data but will not account for the distance between each element of each pair of rows, but for the number of different elements only.
#'
#' @param x A matrix or a data frame.
#' @param byrow By default (\code{byrow=TRUE}), the function computes differences between rows, but it also can do so by columns (\code{byrow=FALSE}).
#'
#' @return
#' An object of class \code{"dist"}.
#'
#' @importFrom stats as.dist
#'
#' @examples
#' set.seed(1)
#' A<-sample(letters[1:10],20,replace=TRUE) # Set a vector of 20 random letters between 'a' and 'j'
#' A
#' B<-replace(A,c(3,7,9),"k") # Replace three values by a 'k', so we can expect 3 differences between A and B
#' length(which(A!=B)) # We're OK
#' C<-A # Set an identical vector to A, so we can expect 0 differences between A and C (and 3 between B and C)
#' length(which(A!=C)) # OK too
#' D<-replace(A,c(1:4,6:8,11:19),letters[11:26]) # Replace 16 values of A by the 16 following letters
#' length(which(A!=D)) # Still OK
#' M<-matrix(c(A,B,C,D),nrow=20,ncol=4)
#' colnames(M)<-c("A","B","C","D")
#' DF<-data.frame(A,B,C,D)
#' t(M)
#' chardist(t(M))
#' chardist(M,byrow=FALSE)
#' chardist(DF,byrow=FALSE)
#'
#' @export

chardist = function(x,byrow=TRUE){
  if(byrow==FALSE){
    x<-t(x)
  }
  D<-matrix(NA, nrow(x), nrow(x))
  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
      if(j<i){
        d = length(which(x[i,]!=x[j,]))
        D[i,j] = d
      }
    }
  }
  if(!is.null(rownames(x))){
    rownames(D)<-rownames(x)
    colnames(D)<-rownames(x)
  }
  return(as.dist(D))
}
