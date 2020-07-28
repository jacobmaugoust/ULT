#' @title Print function for results of multkw functions (multkw, multkw.m and multkw.perm)
#'
#' @description
#' The multivariate Kruskal-Wallis test functions (multkw, multkw.m and multkw.perm) output an object of multkw.output class.
#' It can be saved to an object, but if not, result is automatically printed using this function.
#' If saved, just do print.multkw.output(object) to see its results.
#'
#' @param x An object of class "multkw.output"
#'
#' @method print multkw.output
#' @export

print.multkw.output<-function(x){
  type<-NULL
  if(length(x)==13){type<-"perm"}
  else{
    if(length(x)==10){type<-"NA"}
    else{
      if(length(x)==5){type<-"stock"}
      else{warning("x is not an output of multkw, multkw.m or multkw.perm functions")}
    }
  }
  cat("\n")
  cat(
    if(type=="perm"|type=="NA"){"Extended Multivariate (EM)"}
    else{"Multivariate"},
    "Kruskal-Wallis rank sum test",
    if(type=="perm"){paste0("with ",x$nmc," Monte-Carlo permutations")},
    "\n\n")
  cat("multivariate data:   ",x$y,"\n",sep="")
  cat("           groups:   ",x$group,"\n\n",sep="")
  if(type=="stock"){cat("multivariate Kruskal-Wallis chi-squared = ",x$test.statistic,", df = ",x$df,", p-value = ",x$p.value,sep="")}
  else{
    if(x$pattern.number==1){
      cat("EM Kruskal-Wallis chi-squared = ",x$W2.c,", df = ",x$nu.c,", p-value = ",x$p.multkw.c.chi2,
          if(type=="perm"){paste0(", p-value with ",x$nmc," permutations = ",x$p.multkw.c.perm)},
          "\n\n",sep="")
    }
    else{
      cat("for complete rows only :","\n")
      cat("EM Kruskal-Wallis chi-squared = ",x$W2.c,", df = ",x$nu.c,", p-value = ",x$p.multkw.c.chi2,
          if(type=="perm"){paste0(", p-value with ",x$nmc," permutations = ",x$p.multkw.c.perm)},
          "\n\n",sep="")
      cat("for all rows except those with only missing data :","\n")
      cat("EM Kruskal-Wallis chi-squared = ",x$W2.m,", df = ",x$nu.m,", p-value = ",x$p.multkw.m.chi2,
          if(type=="perm"){paste0(", p-value with ",x$nmc," permutations = ",x$p.multkw.m.perm)},
          "\n\n",sep="")
    }
  }
}
