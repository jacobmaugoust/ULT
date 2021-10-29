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
#' @param simple_multicol_output Logical. The output of the function if you specify a front and a back colors. Set by default to TRUE, meaning the only output will be the "mixed" color. Otherwise, it will return a list of two elements containing (1) the front, back and mixed colors and (2) the mixed transparency.
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
#' transparent.colors(front_color=random_color,front_alpha=random_alpha_value, output="RGB 255 code")
#' # The returning RGB code corresponds to the color in the plot.
#' # For two colors overlapping
#' random_alpha_value<-runif(1,0,1)
#' random_colors<-c()
#' for(i in 1:2){random_colors<-c(random_colors,rgb(red=runif(1,0,1),green=runif(1,0,1),blue=runif(1,0,1),alpha=random_alpha_value))}
#' plot(1:8,1:8,type="n")
#' points(4:5,4:5,cex=20,pch=21,col=NA,bg=random_colors)
#' # The code for all colors
#' transparent.colors(front_color=random_colors,front_alpha=random_alpha_value,output = "RGB 255 code",simple_multicol_output = FALSE)
#' # For several colors successively overlapping
#' random_alpha_value<-runif(1,0,1)
#' random_colors<-c()
#' for(i in 1:5){random_colors<-c(random_colors,rgb(red=runif(1,0,1),green=runif(1,0,1),blue=runif(1,0,1),alpha=random_alpha_value))}
#' dev.new(width=8.958333,height=6.479167,unit="in",noRStudioGD = TRUE)
#' plot(1:1,1:1,type="n")
#' for(i in 1:5){points((1+(i-1)/15),1,cex=c(60-8*i),pch=21,col=NA,bg=random_colors[i])}
#' only_colors<-matrix(nrow=3,ncol=0,NA)
#' superimposed_colors<-matrix(nrow=3,ncol=0,NA)
#' for(i in 1:5){
#'   only_colors<-cbind(only_colors,transparent.colors(front_color = random_colors[i],front_alpha=random_alpha_value,output="RGB 255 code"))
#'   if(i>1){
#'     superimposed_colors<-cbind(superimposed_colors,transparent.colors(front_color=random_colors[i],front_alpha=random_alpha_value,back_color=rgb(t(superimposed_colors[,i-1]/255)),back_alpha = 1,output="RGB 255 code"))
#'   }
#'   else{
#'     superimposed_colors<-cbind(superimposed_colors,only_colors[,1])
#'   }
#' }
#' # Each transparent color (visible on the top right of each circle, on the top of the first one and the right side of the last one)
#' only_colors
#' # Each successive superimposed color (the aggregation of all colors)
#' superimposed_colors
#' # Demonstrating the interest of transparent.colors, especially while dealing with superimposed transparent colors one has to legend
#' cols<-c("black","gray10","gray20","gray30","gray40","gray50","gray60","pink","magenta","red")
#' alpha<-0.2
#' new_cols<-character(length=10)
#' new_cols[10]<-transparent.colors(front_color = cols[10],front_alpha=alpha)
#' for(i in 9:1){
#'   new_cols[i]<-transparent.colors(front_color=cols[i],front_alpha=alpha,back_color=new_cols[i+1],back_alpha=1)
#' }
#' par(mfrow=c(1,2))
#' plot(1:10,1:10,type="n",main="using 'transparent.colors'",axes=FALSE,bty="n",xlab="",ylab="")
#' for (i in 10:1){
#'   polygon(c(1,1:10,10),c(1,seq(from=(5+i/2),to=4,length.out = 10),1),col=new_cols[i],border=NA)
#' }
#' legend("topright",legend=letters[1:10],pch=21,pt.cex=1.5,col=NA,pt.bg=new_cols,text.col="black",bty="n")
#' plot(1:10,1:10,type="n",main="using 'alpha' in the plot, not\nbeing able to correct this in the legend",axes=FALSE,bty="n",xlab="",ylab="")
#' for (i in 10:1){
#'   polygon(c(1,1:10,10),c(1,seq(from=(5+i/2),to=4,length.out = 10),1),col=scales::alpha(cols[i],alpha),border=NA)
#' }
#' legend("topright",legend=letters[1:10],pch=21,pt.cex=1.5,col=NA,pt.bg=cols,text.col="black",bty="n")

#'
#' @importFrom grDevices col2rgb rgb
#'
#' @export

transparent.colors<-function(front_color,back_color=NA,front_alpha,back_alpha=NA,whole_background=NA,output,simple_multicol_output=TRUE){
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
    if(simple_multicol_output==TRUE){
      return(transparent_colors[[3]])
    }
    else{
      return(list(transparent_colors,mixed_transparency))
    }
  }
}
