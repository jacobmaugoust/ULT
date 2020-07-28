#' @title
#' Function to give a vector of 22 contrasting colors
#'
#' @description
#' This function is simply a vector, as rainbow() is for instance. It returns 20 simple and distinct colors found by Sasha Trubetskoy + black and white.
#'
#' @seealso
#' For further details on this set of 20 of the 22 colors, see \url{https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/}
#'
#' @importFrom gtools ask
#'
#' @param black.and.white Logical. Choose ot include or not black and white and theend of the list of colors. By default set to TRUE.
#'
#' @export

contrasting_palette<-function(value,black.and.white=TRUE){
  colors<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  ncol<-22
  if(black.and.white==FALSE){
    ncol<-20
    colors<-colors[1:20]
  }
  if(missing(value)){
    return(colors)
  }
  else{
    if(value>ncol){
      repeat_colors<-ask(paste0("Only ",ncol," colors are available, do you want to choose less than 20 colors (type the number) or do you agree that some colors to repeat (type 'Y', 'YES' or 'yes')?"))
      if(suppressWarnings(is.na(as.numeric(repeat_colors)))==FALSE){
        repeat_colors<-as.numeric(repeat_colors)
        if(repeat_colors<=ncol){
          value<-repeat_colors
        }
        else{
          warning("You did not what I said...")
        }
        return(colors[1:value])
      }
      else{
        iterations<-value%/%ncol
        additions<-value%%ncol
        return(c(rep(colors,iterations),colors[1:additions]))
      }
    }
    else{
      return(colors[1:value])
    }
  }
}
