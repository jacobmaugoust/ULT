#' @title
#' Function to provide the color code of a transparent color
#'
#' @description
#' This function gives the color code (either the name or the RGB code of the color) of a color given its transparency, with or without a specified background.
#' It can also give the color codes of two colors given their transparency, with or without specified background, and of both transparent colors superimposed.
#' This function is particularily useful to plot transparent colors and to use their "real" transparent value to use in a legend.
#'
#' @details
#' The function can return the color code of a single color with a given transparency and a given background.
#' It can also return the color code of two colors with each a given transparency, and a given background. In this case, it can also return the transparency of the two colors superimposed.
#'
#' @param front_color The color (without alpha value) to be transparent over the background. Can be a name, a RGB code, and RGB values (between 0 and 1, 0 and 100 or 0 and 255).
#' @param back_color The color (without alpha value) to be transparent under the front color and over the background. Can be a name, a RGB code, and RGB values (between 0 and 1, 0 and 100 or 0 and 255). Optional.
#' @param front_alpha The transparency (or alpha value) of the front color. Can be between 0 and 1 and 0 and 100.
#' @param back_alpha The transparency (or alpha value) of the back color. Can be between 0 and 1 and 0 and 100. Optional.
#' @param whole_background The background color to be under the front color and potentially the back color. Cannot be transparent. By default, set to white.
#' @param output The output of the color code(s). Can be "color name" (default value), being the color code in R, a "RGB 255 code" (values between 0 and 255) or a "RGB \% code" (values between 0 and 1).
#'
#' @return
#' If there is a single color to be transparent, the color code or name.
#' If there are two colors to be transparent, returns a list: the first element will contain the color codes or names, the second will return the transparency of the two colors superimposed.
#' The output of this function can be used in a plot (or such) function to define a color (see example).
#'
#' @examples
#' # For a single color to be transparent
#' random_alpha_value<-runif(1,0,1)
#' random_color<-rgb(red=runif(1,0,1),green=runif(1,0,1),blue=runif(1,0,1),alpha=random_alpha_value)
#' plot(1:1,type="n")
#' points(1,1,cex=20,pch=21,col=NA,bg=random_color)
#' transparent_colors(front_color=random_color,front_alpha=random_alpha_value, output="RGB 255 code")
#' # The returning RGB code corresponds to the color in the plot.
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @export

transparent_colors<-function(front_color,back_color=NA,front_alpha,back_alpha=NA,whole_background=NA,output){
  if(all(is.na(back_color)==FALSE)){
    if(all(is.na(back_alpha)==TRUE)){
      stop("Please provide an alpha value, ranging from 0 (fully transparent) to 1, 100 (=100%) or 255 (fully opaque), for each color")
    }
  }
  if(all(is.na(back_alpha)==FALSE)){
    if(all(is.na(back_color)==TRUE)){
      stop("Please remove on of the alpha values or provide a second color")
    }
  }
  if(all(is.na(back_color)&is.na(back_alpha))){
    color<-c()
    if(is.matrix(front_color)==TRUE){
      if(dim(as.matrix(front_color))[1]==3
         &dim(as.matrix(front_color))[2]==1){
        if(names(front_color[1,])=="red"
           &names(front_color[2,])=="green"
           &names(front_color[3,])=="blue"){
          if(front_color>1){
            color<-front_color/255
          }
          else{
            if(all(front_color>=0)){
              color<-front_color
            }
          }
        }
        else{
          warning("Front color is not a RGB matrix, the red/green/blue row names are missing as the output of the col2rgb function; try something like col2rgb(rgb(t(front_color))) or name the rows")
        }
      }
      else{
        warning("Front color has not the size of a RGB matrix as the output of the col2rgb function")
      }
    }
    else{
      if(is.character(front_color)==TRUE){
        color<-col2rgb(front_color)/255
      }
      else{
        warning("Front color is not a RGB code color, a HEX code color or a named color")
      }
    }
    alpha<-c()
    {
      if(is.numeric(front_alpha)==FALSE|length(front_alpha)>1){
        warning("Please set a alpha value, ranging from 0 (fully transparent) to 1, 100 (=100%) or 255 (fully opaque)")
      }
      else{
        if(front_alpha>=0&front_alpha<=1){
          alpha<-front_alpha
        }
        else{
          if(front_alpha>=0&front_alpha<=100){
            alpha<-front_alpha/100
          }
          else{
            if(front_alpha>=0&front_alpha<=255){
              alpha<-front_alpha/255
            }
            else{
              warning("Please set alpha value between 0 and 1, 0 and 100 or 0 and 255")
            }
          }
        }
      }
    }
    if(missing(output)){
      output<-"color name"
    }
    if(missing(whole_background)){
      bg<-c(1,1,1)
    }
    {
      if(output=="color name"){
        transparent_color<-rgb(t((1-alpha)*bg+alpha*color))
        return(transparent_color)
      }
      else{
        if(output=="RGB 255 code"){
          transparent_color<-((1-alpha)*bg+alpha*color)*255
          return(transparent_color)
        }
        else{
          if(output=="RGB % code"){
            transparent_color<-(1-alpha)*bg+alpha*color
            return(transparent_color)
          }
        }
      }
    }
  }
  if(all(is.na(back_color)==FALSE&is.na(back_alpha)==FALSE)){
    colors_input<-list(front_color,back_color)
    colors<-vector(mode="list",length=2)
    for (i in 1:2){
    if(is.matrix(colors_input[[i]])==TRUE){
      if(dim(as.matrix(colors_input[[i]]))[1]==3
         &dim(as.matrix(colors_input[[i]]))[2]==1){
        if(names(colors_input[[i]][1,])=="red"
           &names(colors_input[[i]][2,])=="green"
           &names(colors_input[[i]][3,])=="blue"){
          if(any(colors_input[[i]]>1)){
            colors[[i]]<-colors_input[[i]]/255
          }
          else{
            if(all(colors_input[[i]]>=0)){
              colors[[i]]<-colors_input[[i]]
            }
          }
        }
        else{
          if(i==1){
            warning("Front color is not a RGB matrix, the red/green/blue row names are missing as the output of the col2rgb function; try something like col2rgb(rgb(t(front_color))) or name the rows")
          }
          if(i==2){
            warning("Back color is not a RGB matrix, the red/green/blue row names are missing as the output of the col2rgb function; try something like col2rgb(rgb(t(back_color))) or name the rows")
          }
        }
      }
      else{
        if(i==1){
          warning("Front color is not has not the size of a RGB matrix as the output of the col2rgb function")
        }
        if(i==2){
          warning("Back color is not has not the size of a RGB matrix as the output of the col2rgb function")
        }
      }
    }
    else{
      if(is.character(colors_input[[i]])==TRUE){
        colors[[i]]<-col2rgb(colors_input[[i]])/255
      }
      else{
        if(i==1){
          warning("Front color is not a RGB code color, a HEX code color or a named color")
        }
        if(i==2){
          warning("Back color is not a RGB code color, a HEX code color or a named color")
        }
      }
    }
  }
    alpha_input<-c(front_alpha,back_alpha)
    alpha<-c()
    {
    if(is.numeric(alpha_input)==FALSE|length(alpha_input)>2){
      warning("Please set alpha values, ranging from 0 (fully transparent) to 1, 100 (=100%) or 255 (fully opaque)")
    }
    else{
      if(all(alpha_input>=0&alpha_input<=1)){
        alpha<-alpha_input
      }
      else{
        if(all(alpha_input>=0&alpha_input<=100)){
          alpha<-alpha_input/100
        }
        else{
          if(all(alpha_input>=0&alpha_input<=255)){
            alpha<-alpha_input/255
          }
          else{
            warning("Please set alpha values between 0 and 1, 0 and 100 or 0 and 255")
          }
        }
      }
    }
    }
    if(missing(output)){
      output<-"color name"
    }
    if(missing(whole_background)){
      bg<-c(1,1,1)
    }
    transparent_front_color<-(1-alpha[1])*bg+alpha[1]*colors[[1]]
    transparent_back_color<-(1-alpha[2])*bg+alpha[2]*colors[[2]]
    opaque_mixed_color<-(colors[[1]]*alpha[1]+colors[[2]]*alpha[2]*(1-alpha[1]))/(alpha[1]+alpha[2]*(1-alpha[1]))
    mixed_transparency<-alpha[1]+alpha[2]*(1-alpha[1])
    transparent_mixed_color<-(1-mixed_transparency)*bg+mixed_transparency*opaque_mixed_color
    {
      if(output=="color name"){
        transparent_colors<-c(rgb(t(transparent_front_color)),
                              rgb(t(transparent_back_color)),
                              rgb(t(transparent_mixed_color)))
        names(transparent_colors)<-c("front","back","mixed")
      }
      else{
        if(output=="RGB 255 code"){
          transparent_colors<-list(transparent_front_color*255,
                                   transparent_back_color*255,
                                   transparent_mixed_color*255)
          names(transparent_colors)<-c("front","back","mixed")
        }
        else{
          if(output=="RGB % code"){
            transparent_colors<-list(transparent_front_color,
                                     transparent_back_color,
                                     transparent_mixed_color)
            names(transparent_colors)<-c("front","back","mixed")
          }
        }
      }
    }
    return(list(transparent_colors,mixed_transparency))
  }
}
