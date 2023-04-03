#' @title Split character columns of a matrix/dataframe/table (or a list of it)
#'
#' @description This function outputs an inputted matrix/dataframe/table with splitted given columns provided a token for split and new column names.
#'
#' @param x A (list of) matrix/dataframe/table with columns to split.
#' @param split.col Numeric. The column number(s) to split. Row names can also be splitted if \code{split.col=0}. If \code{x} is a list, the elements should be consistent and the column number to split should be equal for each element.
#' @param split.char Character. The token to split values of the given column(s). If several columns are provided, there can either be a single \code{split.char} if all columns are to be splitted with the same token, or as many \code{split.chars} as columns to split with a token for each. Be aware that for several tokens, one may need to use '\\' before (see examples and \link[base]{regex})
#' @param split.names Character. The name of the final splitted columns. For multiple columns to split, there can be either only two names that are re-used for all splitted columns, or twice new column names as columns to split with a pair for each, or no provided names and the final output will be 'oldcolumnname'-1 and 'oldcolumnname'-2
#' @param list.name Character. If \code{x} is a list, the name of the column name for the names of the elements of \code{x}
#'
#' @examples
#' # Create a matrix with five character vectors, two of each (first and third) having each a token between values (a '-' and a '+')
#' set.seed(1)
#' vec1<-paste(sample(letters[7:9],20,TRUE),toupper(sample(letters[7:9],20,TRUE)),sep="-")
#' vec2<-sample(letters[1:5],20,TRUE)
#' vec3<-paste(sample(letters[1:10],20,TRUE),rep("way",20),sep="+")
#' vec4<-toupper(sample(letters[22:26],20,TRUE))
#' vec5<-sample(c("one","two","three"),20,TRUE)
#' x<-cbind(vec1,vec2,vec3,vec4,vec5)
#' x
#' # Split the first and third columns with their tokens, with default final column names
#' split.cols(x,c(1,3),c("-","\\+"))
#' # Do the same with a bit of variation for three lists
#' x<-list()
#' for(i in 1:3){
#'   set.seed(i)
#'   vec1<-paste(sample(letters[(7:9)+i],20,TRUE),toupper(sample(letters[(7:9)+i],20,TRUE)),sep="-")
#'   vec2<-sample(letters[(1:5)*i],20,TRUE)
#'   vec3<-paste(sample(letters[(1:10)+i],20,TRUE),rep(c("way","path","road")[i],20),sep="+")
#'   vec4<-toupper(sample(letters[(22:26)-(3*i)],20,TRUE))
#'   vec5<-sample(c("one","two","three","four","five")[i:(i+2)],20,TRUE)
#'   x<-c(x,list(cbind(vec1,vec2,vec3,vec4,vec5)))
#' }
#' names(x)<-c("firstest","secondtest","thirdtest")
#' # Now split columns 1 and 3 of each element of x and name it "test row"
#' split.cols(x,c(1,3),c("-","\\+"),list.name="test row")
#'
#' @export split.cols

split.cols<-function(x,split.col=NULL,split.char=NULL,split.names=NULL,list.name="variable"){
  if(is.list(x)){
    out<-do.call("rbind",x)
  }
  else{
    out<-x
  }
  if(!is.null(split.col)){
    if(length(split.char)==1){
      split.char<-rep(split.char,length(split.col))
    }
    newcols<-newcolnames<-c()
    newcolorder<-c(1:ncol(out))
    for(i in 1:length(split.col)){
      if(split.col[i]==0){
        split<-strsplit(rownames(out),split.char[i])
      }
      else{
        split<-strsplit(out[,split.col[i]],split.char[i])
      }
      newcols<-cbind(newcols,cbind(unlist(lapply(split,function(x){x[1]})),unlist(lapply(split,function(x){x[2]}))))

      if(!is.null(split.names)){
        if(length(split.names)==2*length(split.col)){
          newcolnames<-c(newcolnames,split.names[(2*(i-1)+1):(2*(i-1)+2)])
        }
        else{
          if(length(split.names)==2){
            newcolnames<-c(newcolnames,split.names)
          }
          else{
            newcolnames<-c(newcolnames,paste0(colnames(out)[split.col[i]],"-",c(1,2)))
          }
        }
      }
      else{
        newcolnames<-c(newcolnames,paste0(colnames(out)[split.col[i]],"-",c(1,2)))
      }

      if(split.col[i]==0){
        newcolorder<-c(0,0,newcolorder)
      }
      else{
        newcolorder<-c(newcolorder[newcolorder<split.col[i]],rep(split.col[i],2),newcolorder[newcolorder>split.col[i]])
      }
    }
    if(any(newcolorder==0)){
      newcolorder[newcolorder==0]<-ncol(out)+(which(split.col==0)-1)*2+c(1,2)
    }
    for(i in 1:length(split.col)){
      if(split.col[i]==0){next}
      else{
        newcolorder[newcolorder==split.col[i]]<-ncol(out)+(i-1)*2+c(1,2)
      }
    }
    out<-cbind(out,newcols)
    colnames(out)<-c(colnames(out)[colnames(out)!=""],newcolnames)
    out<-out[,newcolorder]
  }
  if(is.list(x)){
    out<-cbind(out,rep(names(x),each=nrow(x[[1]])))
    colnames(out)<-c(colnames(out)[colnames(out)!=""],list.name)
    out<-out[,c(ncol(out),(1:(ncol(out)-1)))]
  }
  rownames(out)<-NULL
  return(out)
}
