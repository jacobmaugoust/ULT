#' @title Plot the mapping of a reconstructed trait averaging several methods and add "thermometers" of methods variation
#'
#' @description This function is a wrapper of the \link[ULT]{plot.mapping} and \link[ULT]{thermo.var} functions.
#' Its goal is to take the arithmetic mean of the values of a trait reconstructed for taxa following different methods and/or method parameters.
#' Such arithmetic mean is considered as the "aggregated" method and is displayed by coloring the tree branches.
#' Additional information about the variation between methods can be displayed by "thermometers" (see also \link[ape]{tiplabels} and \link[ape]{nodelabels}) for each taxon, with a color gradient for the value reconstructed for the array of methods.
#'
#' @param tree The phylogenetic tree to put "thermometers" on
#' @param values The data the colors have to follow; can be a data frame or a matrix with taxa as rows and methods as columns.
#' @param type Optional character. If \code{mapping=TRUE}, whether the values represent values for branches (hence coloring the edge with a single color; \code{type="discrete"}) or taxa (hence coloring the edges with a gradient from a taxon to another; \code{type="continuous"})
#' @param cols.args Optional list. A list of arguments for the reference color palette to be passed to the \link[ULT]{scale.palette} function. These arguments are the palette resolution \code{ncol}, the colors to consider \code{col}, the central color \code{middle.col} if there is (otherwise turning this to \code{NA}), a central value in the \code{values} range \code{middle} (if there is, otherwise turning this to \code{NA}), and the \code{values} steps to follow \code{steps} (if there are, otherwise turning this to \code{NA}). Of these, the parameters \code{ncol} and \code{cols} are the most important; the parameters \code{middle.col} and \code{middle} can be left empty, and the parameter \code{span} is estimated as the range of \code{values} if left empty. By default, a "red-yellow-blue" palette of resolution 100 is computed.
#' @param aggr Optional. The reference \code{values} column for the "aggregated" variable (i.e., the variable taking into account all methods by taking the arithmetic mean of all values for each taxon) that represents the center of the color palette (but not necessarily of the "thermometers"!). The \code{values} data can already contain it as being the last column (\code{aggr=TRUE}, the default), as being absent (\code{aggr=FALSE}, computed by the function), as being one of the columns but no the last one (\code{aggr} being the column name or number refering to it in \code{values}).
#' @param order Optional character. To specify the order of the \code{values} to take into account for the "mapping" and the "thermometers" (if \code{thermo=TRUE}). Default is to consider that \code{values} are sorted in the tips/nodes order (\code{order="phylo"}; 1-Ntip rows of \code{values} being for tips 1-N, and so on for the nodes). Values can also be sorted depending on their names (\code{order="names"}; if the tree AND the values have names for tips AND nodes), according to the tree branches construction (\code{order="edge"}; branches construction is available by asking tree$edge, the numbers refering to tips and nodes), or given a custom order (\code{order} being a vector of the names or of the number of all tips/nodes and of same length than the length of \code{values})
#' @param lims Optional. The type or values of limits for the "mapping" and/or the "thermometers" to consider. If a character, limits can encompass the range of plotted values only (\code{lims="local"}; the average of all methods) or of all values (\code{lims="global"}). If \code{thermo=TRUE}, it can be a vector of length 2 specifying a condition for the "mapping" and the "thermometers" (respectively) or of length 1 specifying the same condition for both. If not a character, it can be a numeric vector of length two specifying the numeric limits to consider (for both the "mapping" and the "thermometers" if \code{thermo=TRUE}) or a list containing two vectors of length two specifying the numeric limits to consider for the "mapping" and for the "thermometers".
#' @param thermo Optional logical. Whether to plot "thermometers" to account for values variation or not. Set to \code{TRUE} by default, automatically turned to \code{FALSE} if the values is a single-columned matrix/dataframe.
#' @param plot.mapping.args Optional list. Arguments to be passed to \link[ULT]{plot.mapping} that are not informed from elsewhere. These are the plot title (\code{title}) and other "mapping" arguments (\code{mapping.args})
#' @param thermo.var.args Optional list. Argmuents to be passed to \link[ULT]{thermo.var} that are not informed from elsewhere. These are the resolution for "thermometers" (\code{resolution}), the choice to plot a bar indicating the location of the aggregated value (\code{aggr.bar}) and its color type (\code{aggr.bar.col}), and various graphical aspects of "thermometers" (their border with \code{border}, their size with \code{cex}, their width with \code{width}, their height with \code{height}, and their position relative to taxa with \code{adj}).
#'
#' @import ape
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get a random tree
#' set.seed(10)
#' tree<-rtree(30)
#' # Get a random distribution of values for tips
#' tipvalues<-matrix(ncol=5,nrow=Ntip(tree))
#' set.seed(20)
#' tipvalues[,1]<-fastBM(tree)
#' # Modify a bit these values, with about one half of variation, to simulate different methods
#' for(i in 2:5){
#'   set.seed(i+11)
#'   tipvalues[,i]<-sapply(tipvalues[,1],function(x){x+runif(1,-0.5,0.5)*diff(range(tipvalues[,1]))})
#' }
#' rownames(tipvalues)<-tree$tip.label
#' # Get the ancestral reconstructions for each "method"
#' ancvalues<-matrix(ncol=5,nrow=Nnode(tree))
#' for(i in 1:5){
#'   set.seed(i*2)
#'   ancvalues[,i]<-fastAnc(tree,tipvalues[,i])
#' }
#' # Collate tips and nodes values
#' values<-rbind(tipvalues,ancvalues)
#' # Map average values (all "methods" for each tip) and get an idea of the variation across "methods"
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE)
#' # Get trait values changes between taxa
#' changes<-apply(values,2,function(x){x[tree$edge[,2]]-x[tree$edge[,1]]})
#' # Do as previously for values changes
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge")
#' # Do the same with a finer color resolution
#' cols.args<-list("fun"="scale.palette","ncols"=1000,"cols"=c("blue","yellow","red"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",cols.args=cols.args)
#' # Use only "local" colors with aggregated bars in "thermometers"
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims="local",thermo.var.args=list(aggr.bar=FALSE))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims="local",thermo.var.args=list(aggr.bar=FALSE),cols.args=cols.args)
#' # Same with "global" colors
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims="global",thermo.var.args=list(aggr.bar=FALSE))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims="global",thermo.var.args=list(aggr.bar=FALSE),cols.args=cols.args)
#' # Sort values with custom order
#' order<-sample(1:nrow(values),nrow(values))
#' # Get aggregated values relying on their names
#' tree$node.label<-as.character((1:Nnode(tree))+Ntip(tree)) # Adding names for nodes
#' rownames(values)<-c(tree$tip.label,tree$node.label) #Adding node labels as rownames to values
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,plot.mapping.args=list("title"="reference"))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,order="names",plot.mapping.args=list("title"="relying on names"))
#' # Rely on the given order, still having names for nodes
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,order=order,plot.mapping.args=list("title"="relying on taxa names order"))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,order=order,plot.mapping.args=list("title"="relying on numeric order with node labels"))
#' # Rely on the given order but without node names
#' tree$node.label<-NULL
#' rownames(values)<-c(tree$tip.label,rep("",Nnode(tree)))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,order=order,plot.mapping.args=list("title"="relying on numeric order without node labels"))
#' # Manually add aggregated values a priori and specify it
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,plot.mapping.args=list("title"="reference"))
#' aggr<-apply(values,1,mean)
#' # First with adding aggregated values as the last column
#' values2<-cbind(values,"aggr"=aggr)
#' aggr.trait.map(tree,values2,type="taxa",aggr=TRUE,plot.mapping.args=list("title"="aggr as default (last column)"))
#' # Second, mixing values columns and specify which is the aggregated one by its name or its number
#' values3<-values2[,sample(1:ncol(values2),ncol(values2))]
#' aggr.trait.map(tree,values2,type="taxa",aggr="aggr",plot.mapping.args=list("title"="aggr somewhere, specifying it by its name"))
#' aggr.trait.map(tree,values3,type="taxa",aggr=which(colnames(values3)=="aggr"),plot.mapping.args=list("title"="aggr somewhere, specifying it by its position"))
#'
#' @export

aggr.trait.map<-function(tree,values,type=c("taxa","branch"),cols.args,aggr=TRUE,order=c("phylo","names","edge"),lims=c("local","global"),thermo=TRUE,plot.mapping.args=NULL,thermo.var.args=NULL){
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

  if(is.character(lims)){
    lims<-lapply(lims,function(x){if(x=="local"){range(values[,aggr])}else if(x=="global"){range(values)}})
    map.lims<-lims[[1]]
    if(thermo){
      if(length(lims)==1){
        thermo.lims<-map.lims
      }
      else{
        thermo.lims<-lims[[2]]
      }
    }
  }
  else if(is.numeric(lims)){
    map.lims<-lims
    if(thermo){
      thermo.lims<-lims
    }
  }
  else if(is.list(lims)){
    map.lims<-lims[[1]]
    thermo.lims<-lims[[2]]
  }

  PM.args<-c(list("tree"=tree,
                "values"=values[,aggr],
                "type"=type,
                "lims"=map.lims,
                "order"=order),
             plot.mapping.args[names(plot.mapping.args)%in%c("mapping.args","title")])
  if(!missing(cols.args)){
    PM.args<-c(PM.args,list("cols.args"=cols.args))
  }
  do.call("plot.mapping",PM.args)

  thermo.lims
  if(thermo){
    root.value<-ifelse(type=="branch",FALSE,TRUE)

    if(length(order)==1&&order=="edge"){
      reorder<-tree$edge[,2]
      reorder[reorder>Ntip(tree)]<-reorder[reorder>Ntip(tree)]-1
      values<-values[order(reorder),]
      order<-"phylo"
    }

    TV.args<-c(list("tree"=tree,
                    "values"=values,
                    "thermo.lims"=thermo.lims,
                    "order"=order,
                    "aggr"=aggr,
                    "root.value"=root.value
                    ),
               thermo.var.args[names(thermo.var.args)%in%c("resolution","aggr.bar","aggr.bar.col","border","cex","width","height","adj")])
    if(!missing(cols.args)){
      if("fun"%in%(names(cols.args))){
        cols.args<-cols.args[-which("fun"%in%names(cols.args))]
      }
      TV.args<-c(TV.args,list("cols.args"=cols.args))
    }

    do.call("thermo.var",TV.args)
  }
}