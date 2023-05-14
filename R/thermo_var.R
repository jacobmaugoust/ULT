#' @title Add "thermometers" of trait variation across reconstruction methods for phylogenetic taxa
#'
#' @description This function adds "thermometers" like the \link[ape]{nodelabels} and \link[ape]{tiplabels} functions. However, these functions only add discrete "thermometers", whereas one can want to draw a gardient of values.
#' This is doable by using the \link[ape]{nodelabels} and \link[ape]{tiplabels} functions, but one has to provide a large matrix of values and a large array of colors (with very small steps) to obtain a continuous-like "thermometer".
#' This function aims to simplify this process, only needing tree, values, and options about the palette of colors to use.
#' This is especially useful when "mapping" at once a continuous trait over a tree according to several methods, in order to get an idea of the local (i.e., taxon-level) variation of the reconstructions.
#' By default, the function considers that a previous mapping has been done and that it aims to do "summarize" various methods; this is considered here by calculating the arithmetic mean value of all methods for each taxon (the median could be inadequate if the number of methods is small).
#' By default, no thermometer is reconstructed for taxa with values that are equal for each method (like ancestral reconstructions of tips.
#'
#' @usage thermo.var(tree,values,cols.args,thermo.lims=c("global","local","asym","sym0","symx"),
#'            resolution=100,order=c("phylo","names"),aggr=TRUE,
#'            aggr.bar=TRUE,aggr.bar.col=c("adaptative","blackwhite","black","anycolor"),
#'            root.value=FALSE,border=FALSE,cex=NULL,width=NULL,height=NULL,adj=c(0.5,0.5))
#'
#' @param tree The phylogenetic tree to put "thermometers" on
#' @param values The data the colors have to follow; can be a data frame or a matrix with taxa as rows and methods as columns.
#' @param cols.args Optional list. A list of arguments for the reference color palette to be passed to the \link[ULT]{scale.palette} function. These arguments are the palette resolution \code{ncol}, the colors to consider \code{col}, the central color \code{middle.col} if there is (otherwise turning this to \code{NA}), a central value in the \code{values} range \code{middle} (if there is, otherwise turning this to \code{NA}), and the \code{values} steps to follow \code{steps} (if there are, otherwise turning this to \code{NA}). Of these, the parameters \code{ncol} and \code{cols} are the most important; the parameters \code{middle.col} and \code{middle} can be left empty, and the parameter \code{span} is estimated as the range of \code{values} if left empty. By default, a "blue-yellow-red" palette of resolution 100 is computed.
#' @param thermo.lims Optional. Upper and lower limits of \code{values} to be considered for the global tree-level color palette. Can be a character to set to the range of plotted values only (\code{lims="local"}; the average of all methods) or of all values (\code{lims="global"}) asymmetrically (taking natural values range, \code{"asy"}), symmetrically around zero (taking further value from zero and its opposite, \code{"sym0"}), or around another value (\code{"symx"}, the arithmetic mean by default); hence it can be a vector of length 1 (choosing limits across or within variables), 2 (adding the asymmetric or symmetric choice), or 3 (adding the central value if symmetric to a given value). Otherwise, a numeric of length 2 explicitly specifying values to take as limits.
#' @param resolution Optional numeric. The resolution of the gradient of each "thermometer", i.e., the number of steps (and of successive colors). By default set to 100.
#' @param order Optional character. To specify the order of the \code{values}. Default is to consider that \code{values} are sorted in the tips/nodes order (\code{order="phylo"}; 1-Ntip rows of \code{values} being for tips 1-N, and so on for the nodes). Values can also be sorted depending on their names (\code{order="names"}; if the tree AND the values have names for tips AND nodes) or given a custom order (\code{order} being a vector of the names of all tips/nodes of same length than the number of rows of \code{values})
#' @param aggr Optional. The reference \code{values} column for the "aggregated" variable (i.e., the variable taking into account all methods by taking the arithmetic mean of all values for each taxon) that represents the center of the color palette (but not necessarily of the "thermometers"!). The \code{values} data can already contain it as being the last column (\code{aggr=TRUE}, the default), as being absent (\code{aggr=FALSE}, computed by the function), as being one of the columns but no the last one (\code{aggr} being the column name or number refering to it in \code{values}).
#' @param aggr.bar Optional logical. Whether to plot a horizontal bar in the "thermometer" to highlight the position of the aggregated value on the "thermometer" of each taxon  (the mean not being necessarily in the middle). Default is set to \code{TRUE}
#' @param aggr.bar.col Optional character. If a horizontal bar is drawn for the position of the "aggregated" value, it can be of several colors. By default, it is the RGB inverse (negative) of the aggregated color (\code{aggr.col="adaptative"}), but it can also be the further black or white color to the aggregated color (\code{aggr.col="blackwhite"}), a single color to use for all "thermometers", or a vector of colors to use for all "thermometers" (with one color for each, the vector being recycled if being of different length than the number of rows of \code{values})
#' @param root.value Optional logical. A value for the root is not always available. This can be specified if there is a root value (\code{root.value=TRUE}) or not (\code{root.value=FALSE}, the default)
#' @param border Optional logical. Whether to plot a black border for "thermometer"s if set to \code{TRUE} or not if set to \code{FALSE} (the default)
#' @param cex Optional numeric. The global size ('Character EXpansion'; see the \link[ape]{nodelabels} and \link[ape]{tiplabels} help pages) of the "thermometers"
#' @param width Optional numeric. The width of the "thermometers" (see the \link[ape]{nodelabels} and \link[ape]{tiplabels} help pages)
#' @param height Optional numeric. The height of the "thermometers" (see see the \link[ape]{nodelabels} and \link[ape]{tiplabels} help pages)
#' @param adj Optional numeric. The horizontal and vertical position of the "thermometers" relative to the taxa (see see the \link[ape]{nodelabels} and \link[ape]{tiplabels} help pages)
#'
#' @import ape
#' @importFrom stats setNames
#' @import graphics
#' @importFrom scales rescale
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get a random tree
#' set.seed(2)
#' tree<-rtree(30)
#' # Get a random distribution of values for tips
#' tipvalues<-matrix(ncol=5,nrow=Ntip(tree))
#' set.seed(2)
#' tipvalues[,1]<-fastBM(tree)
#' # Modify a bit these values, with about one half of variation, to simulate different methods
#' for(i in 2:5){
#'   set.seed(i+1)
#'   tipvalues[,i]<-sapply(tipvalues[,1],function(x){x+runif(1,-0.5,0.5)*diff(range(tipvalues[,1]))})
#' }
#' rownames(tipvalues)<-tree$tip.label
#' # Get the ancestral reconstructions for each "method"
#' ancvalues<-matrix(ncol=5,nrow=Nnode(tree))
#' for(i in 1:5){
#'   set.seed(i)
#'   ancvalues[,i]<-fastAnc(tree,tipvalues[,i])
#' }
#' # Collate tips and nodes values
#' values<-rbind(tipvalues,ancvalues)
#' # Get average values for each taxon
#' values<-cbind(values,"aggr"=apply(values,1,mean))
#' # Do a continuous mapping of the "methods" average
#' plot.mapping(tree,values[,6],title="'negative' aggr bars")
#' # Add simple "thermometers", just specifying the tree, the values, and that there is a value for the root
#' thermo.var(tree,values,root.value=TRUE)
#' # Modify the horizontal bars showing the average values so that they are in black or white
#' plot.mapping(tree,values[,6],title="black/white aggr bars")
#' thermo.var(tree,values,root.value=TRUE,aggr.bar.col="blackwhite")
#' # One can notice that the color of the "thermometer" the bar is at is not necessarily the same than the one used for the mapping; this is because the range of average values is not necessarily the same than the range of all values (i.e., including extreme cases smoothed by averaging values)
#' # Re-plot the mapping with a range of colors based on aggregated values only (left) and on all dataset values (righ) and add "thermometers"
#' par(mfrow=c(1,2))
#' plot.mapping(tree,values[,6],lims=range(values[,6]),title="with aggr limits for map")
#' thermo.var(tree,values,root.value=TRUE,aggr.bar.col="blackwhite",border=TRUE)
#' plot.mapping(tree,values[,6],lims=range(values),title="with global limits for map")
#' thermo.var(tree,values,root.value=TRUE,aggr.bar.col="blackwhite",border=TRUE)
#' par(mfrow=c(1,1))
#'
#' # Now play a bit with the options
#' ## Without bar for the average values
#' plot.mapping(tree,values[,6],title="without aggr bar")
#' thermo.var(tree,values,root.value=TRUE,aggr.bar=FALSE)
#' ## With average values-only color palette AND with "thermometers" indexed on that one
#' ## (plus playing with relative horizontal position of "thermometers")
#' plot.mapping(tree,values[,6],title="with thermos whose limits are those of all values (left) and mapping values (right)")
#' thermo.var(tree,values,root.value=TRUE,thermo.lims="global",adj=0.45) # "Thermometers" with all-methods range as limits
#' thermo.var(tree,values,root.value=TRUE,thermo.lims="local",adj=0.5) #"Thermometers" with aggregated method range as limits
#' ## With average values-only color palette for mapping AND "thermometers"
#' plot.mapping(tree,values[,6],title="aggr colors, no aggr bar",lims=range(values[,6]))
#' thermo.var(tree,values,root.value=TRUE,aggr.bar=FALSE,thermo.lims = "local")
#' ## With global values color palette for mapping AND "thermometers"
#' plot.mapping(tree,values[,6],title="global colors, no aggr bar",lims=range(values))
#' thermo.var(tree,values,root.value=TRUE,aggr.bar=FALSE,thermo.lims = "global")
#' ## Modify colors palette
#' plot.mapping(tree,values[,6],cols=list(fun="scale.palette",cols=c("blue","green3","pink")),title="with other colors",lims=range(values))
#' thermo.var(tree,values,root.value=TRUE,cols.args=list(cols=c("blue","green3","pink")))
#' ## Modify thermo resolution to something very extreme
#' plot.mapping(tree,values[,6],title="with low resolution of 'thermometers' (e.g., tips t15 or t7)")
#' thermo.var(tree,values,root.value=TRUE,resolution=3)
#' ## Specify aggr column
#' plot.mapping(tree,values[,6],title="specifying aggr column")
#' thermo.var(tree,values,root.value=TRUE,aggr=6)
#' ## Compute it during function execution
#' plot.mapping(tree,values[,6],title="running aggr during function")
#' thermo.var(tree,values[,1:5],root.value=TRUE,aggr=FALSE)
#' ## Sort values by their names and specify it
#' tree$node.label<-as.character((1:Nnode(tree))+Ntip(tree)) # Add node labels if wanting to rely on names
#' rownames(values)<-c(tree$tip.label,tree$node.label)
#' plot.mapping(tree,values[,6],title="sorting with names")
#' thermo.var(tree,values[sample(1:nrow(values),nrow(values)),],root.value=TRUE,order="names")
#' ## Sort values with an arbitrary order of taxa labels and specify it
#' plot.mapping(tree,values[,6],title="custom label sorting")
#' set.seed(1)
#' custom.order<-sample(1:nrow(values),nrow(values))
#' thermo.var(tree,values[custom.order,],root.value=TRUE,order=c(tree$tip.label,tree$node.label)[custom.order])
#' ## Sort values with an arbitrary order of taxa and specify it
#' plot.mapping(tree,values[,6],title="custom sorting")
#' thermo.var(tree,values[custom.order,],root.value=TRUE,order=custom.order)
#' tree$node.label<-NULL
#' rownames(values)<-c(tree$tip.label,rep("",Nnode(tree)))
#' plot.mapping(tree,values[,6],title="custom sorting without node names")
#' thermo.var(tree,values[custom.order,],root.value=TRUE,order=custom.order)
#' ## Making "thermometers" a bit bigger
#' plot.mapping(tree,values[,6],title="big thermometers")
#' thermo.var(tree,values,root.value=TRUE,cex=2)
#' ## Making "thermometers" a bit wider
#' plot.mapping(tree,values[,6],title="wide thermometers")
#' thermo.var(tree,values,cols.args=list(cols=c("blue","yellow","red")),root.value=TRUE,width=0.15)
#' ## Making "thermometers" a bit higher
#' plot.mapping(tree,values[,6],title="tall thermometers")
#' thermo.var(tree,values,cols.args=list(cols=c("blue","yellow","red")),root.value=TRUE,height=3)
#'
#' @export

thermo.var<-function(tree,values,cols.args,thermo.lims=c("global","local","asym","sym0","symx"),resolution=100,order=c("phylo","names"),aggr=TRUE,aggr.bar=TRUE,aggr.bar.col=c("adaptative","blackwhite","black","anycolor"),root.value=FALSE,border=FALSE,cex=NULL,width=NULL,height=NULL,adj=c(0.5,0.5)){
  if(!nrow(values)%in%c(Ntip(tree),Nnode(tree),(Ntip(tree)+Nnode(tree))-ifelse(root.value,0,1))){
    stop("Please provide values as a list whose length is equal to the number of tips, of nodes, or of both")
  }

  if(length(order)>1&&any(c("phylo","names")%in%order)){
    order<-"phylo"
  }

  if(is.logical(aggr)){
    if(aggr==TRUE){
      aggr<-ncol(values)
    }
    else if(!aggr){
      values<-cbind(values,"Aggregated"=apply(values,1,mean))
      aggr<-ncol(values)
    }
  }
  else if(is.character(aggr)&&as.character(aggr)%in%colnames(values)){
    aggr<-which(colnames(values)==as.character(aggr))
  }
  else if(is.numeric(aggr)&&aggr>ncol(values)){
    stop("Please provided a column number for aggregated values within the number of columns of the values dataset, or any interpretable value for aggr; see help")
  }

  char.to.num.lims<-function(lims,data,i){
    a<-lims[lims%in%c("local","global")][1]
    b<-lims[lims%in%c("asym","sym0","symx")][1]
    x<-lims[!lims%in%c("local","global","asym","sym0","symx")][1]
    if(is.character(x)){x<-as.numeric(x)}

    if(is.na(a)){a<-"local"}
    if(is.na(b)){b<-"asym"}

    if(a=="local"){
      data<-data[,i]
    }
    if(is.na(x)){x<-mean(range(data))}
    if(b!="asym"){
      lims<-c(-1,1)*(max(abs(data-ifelse(b=="symx",x,0))))+ifelse(b=="symx",x,0)
    }
    else{
      lims<-range(data)
    }
    lims
  }
  if(all(is.character(thermo.lims))){
    if("local"%in%thermo.lims){
      lims<-"local"
    }
    else{
      lims<-"global"
    }
    thermo.lims<-do.call("char.to.num.lims",list(thermo.lims,values,aggr))
  }
  else if(is.list(thermo.lims)&&length(thermo.lims)==2){
    lims<-thermo.lims[[which(unlist(lapply(thermo.lims,is.character)))]]
    thermo.lims<-thermo.lims[[which(unlist(lapply(thermo.lims,is.numeric)))]]
  }

  if(missing(cols.args)){
    cols.args<-list("cols"=c("blue","yellow","red"))
  }

  if(length(cols.args)==2&&"hidden"%in%names(cols.args)){
    branch.col.freqs<-cols.args$hidden$branch.col.freqs
    branch.col.freqs.type<-cols.args$hidden$branch.col.freqs.type
    cols.args<-cols.args$cols.args
    if(is.null(cols.args)){
      cols.args<-list("cols"=c("blue","yellow","red"))
    }
  }
  else{
    branch.col.freqs<-branch.col.freqs.type<-NULL
  }

  if(!is.list(cols.args)&!is.null(branch.col.freqs)){
    cols.args<-list("cols"=cols.args)
  }

  if(is.list(cols.args)){
    if(!"ncols"%in%names(cols.args)){cols.args$ncols<-max(1000,resolution*10)}
    if(!"middle.col"%in%names(cols.args)){cols.args$middle.col<-NA}
    if(!"middle"%in%names(cols.args)){cols.args$middle<-NA}
    if(!"span"%in%names(cols.args)){cols.args$span<-thermo.lims}
    if(!is.null(branch.col.freqs)){
      cols<-freq.cols(cols=do.call("scale.palette",cols.args),ncols=cols.args$ncols,freqs=branch.col.freqs,lims=thermo.lims,type=branch.col.freqs.type,values=values[,if("global"%in%lims){c(1:ncol(values))}else{aggr}])
    }
    else{
      cols<-do.call("scale.palette",cols.args)
    }
  }
  else{
    cols<-cols.args
  }

  if(any(aggr.bar.col%in%c("adaptative","blackwhite"))){
    aggr.bar.col<-aggr.bar.col[aggr.bar.col%in%c("adaptative","blackwhite")][1]
  }
  else if(length(aggr.bar.col)!=nrow(values)){
    aggr.bar.col<-rep(aggr.bar.col,length.out=nrow(values))
  }

  taxa<-c(tree$tip.label,if(all(!is.null(tree$node.label))){tree$node.label}else{as.character(c(1:Nnode(tree))+Ntip(tree))})
  check.taxa<-rep(FALSE,length(taxa))

  if(nrow(values)%in%c(Ntip(tree),(Nnode(tree)+Ntip(tree))-ifelse(root.value,0,1))){
    check.taxa[1:Ntip(tree)]<-TRUE
  }
  if(nrow(values)%in%c(Nnode(tree),(Nnode(tree)+Ntip(tree))-ifelse(root.value,0,1))){
    check.taxa[(1:Nnode(tree))+Ntip(tree)]<-TRUE
  }
  if(!root.value){
    check.taxa[(Ntip(tree)+1)]<-FALSE
  }

  parusr <- par("usr")
  if(is.null(cex)){
    cex<-0.75
  }
  if (is.null(width)) {
    width <- cex * (parusr[2] - parusr[1]) / 75
  }
  if (is.null(height)) {
    height <- cex * (parusr[4] - parusr[3]) / 20
  }

  if (length(adj) == 1) {
    adj <- c(adj, 0.5)
  }

  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x<-setNames(lastPP$xx,taxa)[check.taxa]
  y<-setNames(lastPP$yy,taxa)[check.taxa]

  if(!(length(order)==1&&order=="phylo")){
    VM<-FALSE
    if(length(order)==1){
      ref<-rownames(values)
    }
    else if(length(order)!=nrow(values)){
      warn.msg<-"Order provided is not of same length than number of values; considering "
      if(!is.null(rownames(values))){
        ref<-rownames(values)
        warning(paste0(warn.msg,"values row names"))
      }
      else{
        ref<-c(1:nrow(values))
        warning(paste0(warn.msg,"them in same order than tips and nodes"))
      }
    }
    else{
      if(is.character(order)){
        ref<-order
      }
      else if(is.numeric(order)){
        values.match<-order(order)
        VM<-TRUE
      }
    }

    if(!VM){
      if(!is.null(tree$node.label)){
        values.match<-match(c(tree$tip.label,tree$node.label),ref)
      }
      else{
        values.match<-c(match(tree$tip.label,ref[ref%in%tree$tip.label]),which(!ref%in%tree$tip.label))
      }
    }

    values<-values[values.match,,drop=FALSE]
  }

  if(!root.value){
    check.taxa<-check.taxa[-(Ntip(tree)+1)]
  }
  values<-values[which(check.taxa),,drop=FALSE]

  col.values<-apply(values,c(1,2),function(x){
    colnumber<-scales::rescale(x,c(0.5+1e-10,cols.args$ncols+0.5-1e-10),thermo.lims)
    if(colnumber<1){colnumber<-1}
    if(colnumber>cols.args$ncols){colnumber<-cols.args$ncols}
    cols[colnumber]
    })

  for(i in 1:length(x)){
    if(length(unique(values[i,]))==1){next}

    xl <- x[i] - width/2 + adj[1] - 0.5
    xr <- xl + width
    yb <- y[i] - height/2 + adj[2] - 0.5
    yt <- yb + height
    ydiff<-height/resolution

    local.cols<-col.values[i,order(values[i,])]
    rect.scale<-scale.palette(ncols=resolution,cols=c(local.cols[1],local.cols,local.cols[length(local.cols)]),middle.col=col.values[i,aggr],
                              middle=which(names(local.cols)==colnames(col.values)[aggr]),span=c(0.5,length(local.cols)+0.5),steps=c(1:length(local.cols)))

    for(j in 1:length(rect.scale)){
      rect(xl,yb+(j-1)*ydiff,xr,yb+j*ydiff,border=NA,col=rect.scale[j])
    }

    if(aggr.bar){
      barcol<-ifelse(aggr.bar.col=="adaptative",do.call("rgb",c((as.list(apply(col2rgb(local.cols[which(names(local.cols)==colnames(col.values)[aggr])]),1,function(x){255-x}))),list("maxColorValue"=255))),
                     ifelse(aggr.bar.col=="blackwhite",ifelse(mean(col2rgb(col.values[i,aggr]))>0.5*255,"black","white"),aggr.bar.col))
      bary<-yb+ydiff*resolution*which(names(local.cols)==colnames(col.values)[aggr])/(1+length(local.cols))
      segments(xl,bary,xr,bary,col=barcol,lwd=2)
    }

    if(border){
      rect(xl,yb,xr,yt,border="black",col=NA)
    }
  }
}
