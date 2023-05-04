#' @title Plot the mapping of a reconstructed trait averaging several methods and add "thermometers" of methods variation
#'
#' @description This function is a wrapper of the \link[ULT]{plot.mapping} and \link[ULT]{thermo.var} functions.
#' Its goal is to take the arithmetic mean of the values of a trait reconstructed for taxa following different methods and/or method parameters.
#' Such arithmetic mean is considered as the "aggregated" method and is displayed by coloring the tree branches.
#' Additional information about the variation between methods can be displayed by "thermometers" (see also \link[ape]{tiplabels} and \link[ape]{nodelabels}) for each taxon, with a color gradient for the value reconstructed for the array of methods.
#'
#' @usage aggr.trait.map(tree,values,type=c("taxa","branch"),plot=c("methods","aggr"),aggr=TRUE,
#'                order=c("phylo","names","edge"),cols.args,lims=c("local","global","asym","sym0","symx"),
#'                disag=FALSE,return.disag=FALSE,thermo=TRUE,plot.mapping.args=NULL,thermo.var.args=NULL)
#'
#' @param tree The phylogenetic tree to put "thermometers" on.
#' @param values The data to "map" onto the phylogenetic tree and that the colors have to follow; can be a data frame or a matrix with taxa as rows and methods as columns.
#' @param type Optional character. If \code{mapping=TRUE}, whether the values represent values for branches (hence coloring the edge with a single color; \code{type="discrete"}) or taxa (hence coloring the edges with a gradient from a taxon to another; \code{type="continuous"})
#' @param plot Optional character. Whether to plot all methods + "aggregated" values (\code{plot=c("methods","aggr")}, the default) or either one or the other(hence automatically not plotting "thermometers" if \code{plot="methods"}).
#' @param aggr Optional. The reference \code{values} column for the "aggregated" variable (i.e., the variable taking into account all methods by taking the arithmetic mean of all values for each taxon) that represents the center of the color palette (but not necessarily of the "thermometers"!). The \code{values} data can already contain it as being the last column (\code{aggr=TRUE}, the default), as being absent (\code{aggr=FALSE}, computed by the function), as being one of the columns but no the last one (\code{aggr} being the column name or number refering to it in \code{values}).
#' @param order Optional character. To specify the order of the \code{values} to take into account for the "mapping" and the "thermometers" (if \code{thermo=TRUE}). Default is to consider that \code{values} are sorted in the tips/nodes order (\code{order="phylo"}; 1-Ntip rows of \code{values} being for tips 1-N, and so on for the nodes). Values can also be sorted depending on their names (\code{order="names"}; if the tree AND the values have names for tips AND nodes), according to the tree branches construction (\code{order="edge"}; branches construction is available by asking tree$edge, the numbers refering to tips and nodes), or given a custom order (\code{order} being a vector of the names or of the number of all tips/nodes and of same length than the length of \code{values})
#' @param cols.args Optional list. A list of arguments for the reference color palette to be passed to the \link[ULT]{scale.palette} function. These arguments are the palette resolution \code{ncol}, the colors to consider \code{col}, the central color \code{middle.col} if there is (otherwise turning this to \code{NA}), a central value in the \code{values} range \code{middle} (if there is, otherwise turning this to \code{NA}), and the \code{values} steps to follow \code{steps} (if there are, otherwise turning this to \code{NA}). Of these, the parameters \code{ncol} and \code{cols} are the most important; the parameters \code{middle.col} and \code{middle} can be left empty, and the parameter \code{span} is estimated as the range of \code{values} if left empty. By default, a "red-yellow-blue" palette of resolution 100 is computed. Please not that if a lot of colors are provided, if \code{type="branch"}, and if \code{branch.col.freqs} and optionally \code{branch.col.freqs.type} are provided in \code{plot.mapping.args}, the function may encounter a bug: some points would not have an attributed color because of the too narrow value steps between each color.
#' @param lims Optional. The type or values of limits for the "mapping" and/or the "thermometers" to consider. If a character, limits can encompass the range of plotted values only (\code{lims="local"}; the average of all methods) or of all values (\code{lims="global"}) and they can be asymmetric (taking natural values range, \code{"asy"}), symmetric around zero (taking further value from zero and its opposite, \code{"sym0"}), or around another value (\code{"symx"}, the arithmetic mean by default); hence it can be a vector of length 1 (choosing limits across or within variables), 2 (adding the asymmetric or symmetric choice), or 3 (adding the central value if symmetric to a given value) that is applicable for all desired "mappings" and for "thermometers. If a numeric, in the case of a continuous trait, the value limits to consider (two values: the inferior and superior bounds), hence a numeric vector of length two. In both cases, it can also be a list of such vectors (character or numeric) of length two if \code{thermo=TRUE} and desiring different limits for "mappings" and "thermometers", or a list of same length than the number of desired plots (depending on what is passed to the \code{plot} argument, with therefore different limits for each plot) and whose elements are either a single vector (specifying limits for "mappings" and also for "thermometers" if applicable) or a list of two vectors (specifying limits for both "mappings" and "thermometers").
#' @param disag Optional. If \code{type="branch"}, whether to consider "disagreements" between branch values; this is especially useful while considering signs or discrete traits. By default set to \code{FALSE}, can be set to \code{TRUE} (hence coloring disagreeing branches in black) or to any color (hence coloring disagreeing branches in the given color).
#' @param return.disag If \code{disag=TRUE}, whether to return the list of disagreeing taxa.
#' @param thermo Optional logical. Whether to plot "thermometers" to account for values variation or not. Set to \code{TRUE} by default, automatically turned to \code{FALSE} if the values is a single-columned matrix/dataframe.
#' @param plot.mapping.args Optional list. Arguments to be passed to \link[ULT]{plot.mapping} that are not informed from elsewhere. These are the plot title (\code{title}) and other "mapping" arguments (\code{mapping.args}). Can be a list of lists if \code{plot="methods"} or \code{plot=c("methods","aggr")}, then of same length than the desired number of plots.
#' @param thermo.var.args Optional list. Argmuents to be passed to \link[ULT]{thermo.var} that are not informed from elsewhere. These are the resolution for "thermometers" (\code{resolution}), the choice to plot a bar indicating the location of the aggregated value (\code{aggr.bar}) and its color type (\code{aggr.bar.col}), and various graphical aspects of "thermometers" (their border with \code{border}, their size with \code{cex}, their width with \code{width}, their height with \code{height}, and their position relative to taxa with \code{adj}).
#'
#' @import ape
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get a random tree
#' set.seed(10)
#' tree<-rtree(50)
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
#' # Map all "methods" + average values (all "methods" for each tip) and get an idea of the variation across "methods"
#' layout(get.grid(ncol(values)+1))
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE)
#' # All plots have their own color range according to their values, but there are obviously some discrepancies, let's use the total data range instead
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims="global")
#' # The range of values being not especially controlled but having values above and below zero, let's have the first condition ("local" color ranges) with now a symmetric color scale
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims="sym0")
#' # Let's now have a symmetrical color range indexed on all values
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims=c("global","sym0"))
#' # The data are not necessarily distributed around zero, let's have a symmetric color range but around the "real" middle of the values range, to have closer fidelity of the colors to the values
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims=c("global","symx",mean(range(values))))
#' # Let's repeat this but this time using "local" color ranges for "mappings" (to have maximal palette range for each), still keeping a symmetric color palette (using thus "symx") centered on "local" data arithmetic mean (not providing the "x" that will be automatically calculated)
#' # and "global" color ranges for "thermometers" (to check the "aggregated" values position, with necessarily the colors of the "thermometers" not matching that of the "aggregated mapping")
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,lims=list(c("local","symx"),c("global","symx",mean(range(values)))))
#'
#' # Get trait values changes between taxa
#' changes<-apply(values,2,function(x){x[tree$edge[,2]]-x[tree$edge[,1]]})
#' # Do as previously for values changes
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge")
#' # This time, we can now that changes are negative or positive, so with a "center" around zero; let's have a symmetric color range around zero
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims="sym0")
#' # Due to a very high extreme we lost a lot of finer resolution, so let's do the same but with
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"))
#' # As it stands, it converts changes values in three colors, hence separating distribution (governed by limits) in three ranges of equal width.
#' # However, we could be interested in having an uneven repartition since most branches are in the second condition (color yellow), i.e., close to a zero value.
#' # Let's have the same as previously but with colors distributed so to have three even classes of colors.
#' # To do so, we need to indicate that the colors do not represent equal thirds between data limits but a custom width (representing, here, one third of values distribution)
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"),plot.mapping.args = list("branch.col.freqs"="equal","branch.col.freqs.type"="proportion"))
#' # Do the same with a finer color resolution
#' cols.args<-list("fun"="scale.palette","ncols"=1000,"cols"=c("blue","yellow","red"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",cols.args=cols.args,lims=c("local","sym0"),plot.mapping.args = list("branch.col.freqs"="equal","branch.col.freqs.type"="proportion"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",cols.args=cols.args,lims=c("global","sym0"),plot.mapping.args = list("branch.col.freqs"="equal","branch.col.freqs.type"="proportion"))
#' # Now do the same with the values in the 5% range around zero being in yellow
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("local","sym0"),plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="width"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"),plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="width"))
#' # Now do the same with the values in the 5% range around zero being in yellow before doing the color gradient (with therefore more colors, the values in the 5% range around zero being the yellowests)
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("local","sym0"),cols.args=cols.args,plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="width"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"),cols.args=cols.args,plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="width"))
#' # Now do the same with the 5% central values being in yellow
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("local","sym0"),plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="proportion"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"),plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="proportion"))
#' # Now do the same with the 5% central values being in yellow before doing the color gradient (with therefore more colors, the 5% central values being the yellowests)
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("local","sym0"),cols.args=cols.args,plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="proportion"))
#' aggr.trait.map(tree,changes,type="branch",aggr=FALSE,order="edge",lims=c("global","sym0"),cols.args=cols.args,plot.mapping.args = list("branch.col.freqs"=c(0.475,0.05,0.475),"branch.col.freqs.type"="proportion"))
#'
#' # Get signs of changes
#' signs<-apply(changes,c(1,2),sign)
#' aggr.trait.map(tree,signs,type="branch",aggr=FALSE,order="edge",disag=TRUE)
#' # Now signs are too much present, let's get some tolerance and set as zero value the ones close to zero
#' range(changes) # Let's take -0.5 and 0.5 as bounds
#' signs<-apply(changes,c(1,2),function(x){if(x<0.5&x>(-0.5)){x<-0}else{sign(x)}})
#' aggr.trait.map(tree,signs,type="branch",aggr=FALSE,order="edge",disag=TRUE)
#' # Let's do same with -1 and 1 as bounds, and with disagreeing branches in light gray, to first see agreeing branches
#' signs<-apply(changes,c(1,2),function(x){if(x<1&x>(-1)){x<-0}else{sign(x)}})
#' aggr.trait.map(tree,signs,type="branch",aggr=FALSE,order="edge",disag="lightgray")
#' # Now let's get the disagreeing branches
#' disag<-aggr.trait.map(tree,signs,type="branch",aggr=FALSE,order="edge",disag="lightgray",return.disag=TRUE)
#' disag[which(disag)]
#'
#' # Sort values with custom order
#' order<-sample(1:nrow(values),nrow(values))
#' # Get aggregated values relying on their names
#' tree$node.label<-as.character((1:Nnode(tree))+Ntip(tree)) # Adding names for nodes
#' rownames(values)<-c(tree$tip.label,tree$node.label) #Adding node labels as rownames to values
#' par(mfrow=c(1,1))
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,plot="aggr",plot.mapping.args=list("title"="reference"))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,plot="aggr",order="names",plot.mapping.args=list("title"="relying on names"))
#' # Rely on the given order, still having names for nodes
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,plot="aggr",order=order,plot.mapping.args=list("title"="relying on taxa names order"))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,plot="aggr",order=order,plot.mapping.args=list("title"="relying on numeric order with node labels"))
#' # Rely on the given order but without node names
#' tree$node.label<-NULL
#' rownames(values)<-c(tree$tip.label,rep("",Nnode(tree)))
#' aggr.trait.map(tree,values[order,],type="taxa",aggr=FALSE,plot="aggr",order=order,plot.mapping.args=list("title"="relying on numeric order without node labels"))
#' # Manually add aggregated values a priori and specify it
#' aggr.trait.map(tree,values,type="taxa",aggr=FALSE,plot="aggr",plot.mapping.args=list("title"="reference"))
#' aggr<-apply(values,1,mean)
#' # First with adding aggregated values as the last column
#' values2<-cbind(values,"aggr"=aggr)
#' aggr.trait.map(tree,values2,type="taxa",aggr=TRUE,plot="aggr",plot.mapping.args=list("title"="aggr as default (last column)"))
#' # Second, mixing values columns and specify which is the aggregated one by its name or its number
#' values3<-values2[,sample(1:ncol(values2),ncol(values2))]
#' aggr.trait.map(tree,values2,type="taxa",aggr="aggr",plot="aggr",plot.mapping.args=list("title"="aggr somewhere, specifying it by its name"))
#' aggr.trait.map(tree,values3,type="taxa",aggr=which(colnames(values3)=="aggr"),plot="aggr",plot.mapping.args=list("title"="aggr somewhere, specifying it by its position"))
#'
#' @export

aggr.trait.map<-function(tree,values,type=c("taxa","branch"),plot=c("methods","aggr"),aggr=TRUE,order=c("phylo","names","edge"),cols.args,lims=c("local","global","asym","sym0","symx"),disag=FALSE,return.disag=FALSE,thermo=TRUE,plot.mapping.args=NULL,thermo.var.args=NULL){
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

  branch<-any("branch"%in%type)

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

  if(!"aggr"%in%plot){
    nplots<-ncol(values)-1
    plot.order<-c(1:nplots)
    thermo<-FALSE
    thermo.var.args<-NULL
  }
  else if(!"methods"%in%plot){
    nplots<-1
    plot.order<-aggr
  }
  else{
    nplots<-ncol(values)
    plot.order<-c(c(1:nplots)[-aggr],aggr)
  }

  if(is.list(lims)){
    if(length(lims)==nplots){
      if(is.list(lims[[1]])){
        map.lims<-lapply(lims,function(x){x[1]})
        thermo.lims<-lapply(lims,function(x){x[2]})
      }
      else{
        map.lims<-lims
        if(thermo){
          thermo.lims<-lims
        }
      }
    }
    else if(length(lims)==2){
      map.lims<-rep(lims[1],nplots)
      thermo.lims<-rep(lims[2],nplots)
    }
  }
  else if(is.vector(lims)){
    map.lims<-rep(list(lims),nplots)
    if(thermo){
      thermo.lims<-rep(list(lims),nplots)
    }
  }

  if(!is.null(plot.mapping.args)&&(length(plot.mapping.args)!=nplots|!is.null(names(plot.mapping.args)))){
    plot.mapping.args<-rep(list(plot.mapping.args),nplots)
  }

  if(type=="branch"&disag!=FALSE&"aggr"%in%plot){
    if(is.character(disag)){
      disag.col<-disag
      disag<-TRUE
    }
    else{
      disag.col<-"black"
    }
  }

  for(i in 1:length(plot.order)){
    if(all(is.character(map.lims[[i]]))){
      temp.map.lims<-do.call("char.to.num.lims",list(map.lims[[i]],values,plot.order[i]))
    }
    else{
      temp.map.lims<-map.lims[[i]]
    }
    PM.args<-c(list("tree"=tree,
                    "values"=values[,plot.order[i]],
                    "type"=type,
                    "lims"=temp.map.lims,
                    "order"=order),
               plot.mapping.args[[i]][names(plot.mapping.args[[i]])%in%c("title","branch.col.freqs","branch.col.freqs.type","mapping.args")])
    if(!missing(cols.args)){
      PM.args<-c(PM.args,list("cols.args"=cols.args))
    }

    if(branch&&"branch.col.freqs"%in%names(PM.args)&&("global"%in%lims|temp.map.lims[1]>min(PM.args$values)&temp.map.lims[2]<max(PM.args$values))){
      PM.args$branch.col.freqs<-list("branch.col.freqs"=PM.args$branch.col.freqs,
                                     "global.values"=values[values>=temp.map.lims[1]&values<=temp.map.lims[2]])
    }

    do.call("plot.mapping",PM.args)

    if(type=="branch"&disag&"aggr"%in%plot&plot.order[i]==aggr){
      disag.values<-apply(values[,-aggr,drop=FALSE],1,function(x){length(unique(sign(x)))>1})
      lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      last_disag<-which(disag.values)
      for(d in 1:length(last_disag)){
        temp_x<-lastPP$xx[tree$edge[last_disag[d],]]
        temp_y<-lastPP$yy[tree$edge[last_disag[d],]]
        segments(temp_x[1],temp_y[1],temp_x[1],temp_y[2],col=disag.col,lwd=3)
        segments(temp_x[1],temp_y[2],temp_x[2],temp_y[2],col=disag.col,lwd=3)
      }
      if(return.disag){
        return(disag.values)
      }
    }


    if(thermo&&plot.order[i]==aggr){
      root.value<-ifelse(type=="branch",FALSE,TRUE)

      if(length(order)==1&&order=="edge"){
        reorder<-tree$edge[,2]
        reorder[reorder>Ntip(tree)]<-reorder[reorder>Ntip(tree)]-1
        values<-values[order(reorder),]
        order<-"phylo"
      }

      TV.args<-c(list("tree"=tree,
                      "values"=values,
                      "thermo.lims"=thermo.lims[[i]],
                      "order"=order,
                      "aggr"=aggr,
                      "root.value"=root.value
                      ),
                 thermo.var.args[names(thermo.var.args)%in%c("resolution","aggr.bar","aggr.bar.col","border","cex","width","height","adj")])

      if(!missing(cols.args)){
        if("fun"%in%(names(cols.args))){
          cols.args<-cols.args[-which(names(cols.args)=="fun")]
        }
        TV.args<-c(TV.args,list("cols.args"=cols.args))
      }
      if(branch&&"branch.col.freqs"%in%names(PM.args)){
        if("global"%in%lims|temp.map.lims[1]>min(PM.args$values)&temp.map.lims[2]<max(PM.args$values)){
          hidden<-PM.args$branch.col.freqs[names(PM.args$branch.col.freqs)=="branch.col.freqs"]
        }
        else{
          hidden<-PM.args[names(PM.args)=="branch.col.freqs"]
        }
        if("branch.col.freqs.type"%in%names(PM.args)){
          hidden<-c(hidden,PM.args[names(PM.args)=="branch.col.freqs.type"])
        }
        TV.args$cols.args<-c("cols.args"=if("cols.args"%in%names(TV.args)){list(TV.args$cols.args)}else{list(NULL)},"hidden"=list(hidden))
      }
      do.call("thermo.var",TV.args)
    }
  }
}
