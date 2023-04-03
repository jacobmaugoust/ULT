#' @title Invert dimensions of a list of lists (with 2 or more successive listings)
#'
#' @description Inverts the successive listings of a list of lists
#'
#' @importFrom stats setNames
#'
#' @param lists The list of lists (with 2 or more successive listings)
#' @param new.order A numeric vector specifying the re-ordering of all sub-lists
#'
#' @examples
#' # Get a list with various type of elements and three sublevels
#' set.seed(1)
#' list<-lapply(c(1:5),function(x){lapply(c(1:4),function(y){lapply(c(1:3),function(z){
#'   n<-sample(c(1:5),1)
#'   out<-NULL
#'   if(n==1){out<-letters[1:sample(c(2:26),1)]}
#'   else if(n==2){out<-runif(sample(c(10:100),1))}
#'   else if(n==3){out<-sample(c(TRUE,FALSE),50,replace=TRUE)}
#'   else if(n==4){out<-matrix("a",ncol=4,nrow=sample(c(6:12),1))}
#'   else if(n==5){out<-"hey"}
#'   out
#' })})})
#' list
#' length(list)
#' length(list[[1]])
#' length(list[[1]][[1]])
#' # Without names to the list
#' invert.dim.lists(list,c(3,2,1))
#' # Set names for each elements: x for first order elements, y for second order ones, z for third order ones (so that we can check whether the function performs well)
#' names(list)<-paste("x",1:length(list),sep="")
#' for(i in 1:length(list)){
#'   names(list[[i]])<-paste("y",1:length(list[[i]]),sep="")
#'   for(j in 1:length(list[[i]])){
#'     names(list[[i]][[j]])<-paste("z",1:length(list[[i]][[j]]),sep="")
#'   }
#' }
#' list
#' invert.dim.lists(list,c(3,2,1))
#' @export

invert.dim.lists<-function(lists,new.order){
  test<-lists
  nL<-0
  while(is.list(test)){
    nL<-nL+1
    test<-test[[1]]
  }
  if(any(new.order>nL)){
    new.order<-new.order[new.order<=nL]
  }
  if(nL>length(new.order)){
    new.order<-c(new.order,c(1:nL)[!c(1:nL)%in%new.order])
  }

  names<-list()
  for(i in 1:nL){
    temp_names<-names(eval(parse(text=paste0("lists",ifelse(new.order[i]==1,"",paste0("[[",rep(1,new.order[i]-1),"]]",collapse=""))))))
    if(is.null(temp_names)){
      names[[i]]<-NA
    }
    else{
      names[[i]]<-temp_names
    }
  }

  res<-NA
  for(i in nL:1){
    res<-setNames(list(rep(list(res),length(eval(parse(text=paste0("lists",ifelse(new.order[i]==1,"",paste0("[[",rep(1,new.order[i]-1),"]]",collapse=""))))))))[[1]],if(!all(is.na(names[[i]]))){names[[i]]})
  }

  expr<-""
  for(i in c(1:nL)){
    expr<-paste0(expr,"for (",letters[i]," in 1:length(lists",ifelse(new.order[i]==1,")",paste0(paste0("[[",rep(1,new.order[i]-1),"]]",collapse=""),")")),"){")
  }
  expr<-paste0(expr,"res",paste0("[[",letters[1:nL],"]]",collapse=""),"<-","lists",paste0("[[",letters[1:nL][order(new.order)],"]]",collapse=""),paste0(rep("}",nL),collapse=""))

  eval(parse(text=expr))

  res
}
