#' @title "Map" a trait: plot its distribution over a phylogenetic tree
#'
#' @description This function plots continuous and discrete trait onto a phylogenetic tree, with values being either on taxa (tips/nodes) or on branches.
#' It is largely a wrapper of the function \link[phytools]{contMap} for the continous "mapping" and of the function \link[ape]{plot.phylo} for the discrete "mapping".
#'
#' @usage plot.mapping(tree,values,type=c("taxa","branch"),cols.args,title=NULL,lims=NULL,
#'              order=c("phylo","names","edge"),mapping.args=NULL)
#'
#' @param tree The phylogenetic tree to "map" the trait on
#' @param values The (univariate) data to "map" onto the phylogenetic tree
#' @param type The unit the values to plot do represent: either these values stand for branches (\code{type="branch"}), or they stand for taxa (\code{type="taxa"}).
#' @param cols.args The colors to use to represent the values. This can be a character vector specifying all colors to use: if \code{type="branch"}, there should be as many colors as branches, and if \code{type="taxa"}, there should be as many colors as is the resolution of the continuous plot by \code{contMap} (1000 by default) plus 1. This can also be a list of arguments to be used to attribute given colors to be passed to \link[ULT]{discrete.palette} or \link[ULT]{scale.palette} (depending if the trait is discrete or continuous, respectively); in any case, the list should at least include a \code{"cols"} vector of desired colors and a \code{"fun"} vector specifying the function to use ("discrete.palette" or "scale.palette"). If \code{type="branch"}, a simple vector of colors can be provided, and the \code{values} are splitted in as many case as there are colors minus one (the latest color being for the highest values), each range being represented by a color (i.e., potentially different values with same color); it is therefore strongly advised to provide a color palette larger than the number of \code{values} in that special case. By default, it is set to a \code{list(fun="scale.palette",cols="blue","yellow","red")}.
#' @param title Optional character. To provide a title to the plot.
#' @param lims Optional numeric. In the case of a continuous trait, the value limits to consider (two values: the inferior and superior bounds).
#' @param order Optional character. To specify the order of the \code{values} to take into account for the "mapping". Default is to consider that \code{values} are sorted in the tips/nodes order (\code{order="phylo"}; 1-Ntip rows of \code{values} being for tips 1-N, and so on for the nodes). Values can also be sorted depending on their names (\code{order="names"}; if the tree AND the values have names for tips AND nodes), according to the tree branches construction (\code{order="edge"}; branches construction is available by asking tree$edge, the numbers refering to tips and nodes), or given a custom order (\code{order} being a vector of the names or of the number of all tips/nodes and of same length than the length of \code{values})
#' @param mapping.args Optional list. List of arguments to take into account for the "mapping". Depending on the \code{type} of trait, these arguments are either to be passed to \link[ape]{plot.phylo} (for discrete trait) or to \link[phytools]{contMap} (for continuous trait)
#'
#' @import ape
#' @importFrom stats setNames
#' @importFrom graphics par
#' @importFrom scales rescale
#'
#' @examples
#' require(ape)
#' require(phytools)
#' # Get a random tree
#' set.seed(1)
#' tree<-rtree(30)
#'
#' # Work with a continuous trait on all taxa
#' ## Get a random continuous trait
#' set.seed(2)
#' tips<-fastBM(tree)
#' ## Get ancestral states of it
#' ancs<-ace(tips,tree,"continuous","REML")$ace
#' ## Collate both
#' x<-c(tips,ancs) # As it is, it is a vector with tip values then node values, so sorted "phylogenetically", according to tips and nodes order
#' ## Plot the trait using a gradual color palette with blue, yellow, and red (the default)
#' plot.mapping(tree,x,"taxa",cols.args=list(fun="scale.palette",cols=c("blue","yellow","red")))
#' ## Plot the trait using a gradual color palette with blue, green, and pink
#' plot.mapping(tree,x,"taxa",cols.args=list(fun="scale.palette",cols=c("blue","green3","pink")))
#' ## Plot the trait with a title
#' plot.mapping(tree,x,"taxa",title="mapping of continuous x")
#' ## Plot the trait with a legend bar
#' plot.mapping(tree,x,"taxa",mapping.args=list(legend=TRUE))
#' ## Plot the trait with different lims, setting the start at 0 (so that all values below zero should be blue)
#' plot.mapping(tree,x,"taxa",lims=c(0,max(x)),mapping.args=list(legend=TRUE))
#' ## Sort the values in an odd way to rely on names and plot the mapping
#' plot.mapping(tree,x,"taxa",title="reference") # The reference mapping
#' tree$node.label<-as.character((1:Nnode(tree))+Ntip(tree)) # Add node labels so that they match the names of x
#' plot.mapping(tree,x[sample(c(1:length(x)),length(x))],"taxa",order="names",title="with names") # With alphabetically sorted values
#' ## Sort the values in an odd but given way
#' custom.sort<-sample(1:length(x),length(x))
#' plot.mapping(tree,x[custom.sort],"taxa",order=names(x)[custom.sort],title="custom sorting with names") # With custom sorting on names (implies that both the three and the values have tip and node labels)
#' plot.mapping(tree,x[custom.sort],"taxa",order=custom.sort,title="custom sorting") # With custom sorting
#' tree$node.label<-NULL
#' plot.mapping(tree,x[custom.sort],"taxa",order=custom.sort,title="custom sorting with no names") # With custom sorting not being able to rely at all on names since tree has no node labels anymore
#'
#' # Work with a discrete trait on all taxa
#' ## Get a random discrete trait
#' set.seed(3)
#' tips<-round(fastBM(tree),0)
#' ## Get ancestral states of it
#' ancs<-apply(ace(tips,tree,"discrete")$lik.anc,1,function(x){sort(unique(tips))[which.max(x)]})
#' ## Collate both
#' x<-c(tips,ancs) # As previously, it is sorted "phylogenetically" as it is now
#' ## Plot the trait
#' plot.mapping(tree,x,"taxa",cols.args=list(fun="discrete.palette",cols=contrasting.palette(length(unique(x))),ncols=1001))
#' legend("topright",legend=sort(unique(x)),lwd=2,col=contrasting.palette(length(unique(x))),bty="n")
#' ## To check color distribution with values
#' x
#' nodelabels()
#'
#' # Work with a trait on branches
#' ## Get a random continuous trait
#' set.seed(4)
#' tips<-fastBM(tree)
#' ## Get ancestral states of it
#' ancs<-ace(tips,tree,"continuous","REML")$ace
#' ## Collate both
#' trait<-c(tips,ancs)
#' ## Get a gradient of colors according to these values
#' cols.args<-scale.palette(ncols=1001,cols=c("blue","yellow","red"),middle=NA,middle.col=NA,span=range(trait))
#' ## Plot it
#' plot.mapping(tree,trait,cols.args=cols.args,mapping.args=list(legend=TRUE))
#' ## Get trait value changes over each branch (in the order of the tree$edge element)
#' x<-apply(tree$edge,1,function(x){trait[x[2]]-trait[x[1]]})
#' ## Plot the trait value changes with the same color gradient than for the trait values themselves
#' plot.mapping(tree,x,"branch",cols.args=cols.args,order="edge")
#' ## Plot the trait value changes by simply splitting them in three colors (the lower half of the changes being in blue, the second in yellow, and the highest value in red)
#' plot.mapping(tree,x,"branch",cols.args=c("blue","yellow","red"),order="edge")
#' ## Changing the order of the trait values
#' plot.mapping(tree,x,"branch",cols.args=c("blue","yellow","red"),order="edge",title="reference")
#' ### Sorting according to phylogeny (order of tips and nodes)
#' plot.mapping(tree,x[order(tree$edge[,2])],"branch",cols.args=c("blue","yellow","red"),order="phylo",title="phylo sorting")
#' ### Sorting according to taxa labels
#' tree$node.label<-as.character((1:Nnode(tree))+Ntip(tree)) # Assigning node labels
#' plot.mapping(tree,setNames(x,c(tree$tip.label,tree$node.label)[tree$edge[,2]]),"branch",cols.args=c("blue","yellow","red"),order="names",title="with names, no sorting")
#' custom.order<-sample(1:length(x),length(x))
#' plot.mapping(tree,setNames(x,c(tree$tip.label,tree$node.label)[tree$edge[,2]])[custom.order],"branch",cols.args=c("blue","yellow","red"),order="names",title="with names, custom sorting 1")
#' plot.mapping(tree,setNames(x,c(tree$tip.label,tree$node.label)[tree$edge[,2]])[custom.order],"branch",cols.args=c("blue","yellow","red"),order=c(tree$tip.label,tree$node.label)[tree$edge[,2]][custom.order],title="with names, custom sorting 2")
#' tree$node.label<-NULL
#' plot.mapping(tree,x[custom.order],"branch",cols.args=c("blue","yellow","red"),order=custom.order,title="no names, custom sorting")
#'
#' @export plot.mapping

plot.mapping<-function(tree,values,type=c("taxa","branch"),cols.args,title=NULL,lims=NULL,order=c("phylo","names","edge"),mapping.args=NULL){
  if(length(order)>1&&length(order)!=length(values)){
    order<-order[order%in%c("phylo","names","edge")][1]
  }
  if(length(type)>1){
    type<-type[c("taxa","branch")%in%type][1]
  }
  branch<-any("branch"%in%type)

  if(missing(cols.args)){
    cols.args<-list(fun="scale.palette",cols=c("blue","yellow","red"))
  }

  if(is.null(lims)){
    lims<-range(values)
  }

  PP.args<-CM.args<-PCM.args<-NULL

  if(!is.null(mapping.args)){
    plotphylo.args<-c("use.edge.length",
                      "node.pos", "show.tip.label",
                      "show.node.label", "edge.width", "edge.lty", "node.color",
                      "node.width", "node.lty", "font", "cex",
                      "adj", "srt", "no.margin", "root.edge", "label.offset", "underscore",
                      "x.lim", "y.lim", "lab4ut", "tip.color",
                      "rotate.tree", "open.angle", "node.depth", "align.tip.label")
    contmap.args<-c("res", "fsize", "ftype", "lwd", "legend", "outline", "sig", "type")
    plotcontmap.args<-c("legend","fsize","ftype","outline","lwd","sig","type","mar","offset","xlim","ylim","hold","leg.txt")
    PP.args<-mapping.args[names(mapping.args)%in%plotphylo.args]
    CM.args<-mapping.args[names(mapping.args)%in%contmap.args]
    PCM.args<-mapping.args[names(mapping.args)%in%plotcontmap.args]
  }

  if(is.list(cols.args)){
    if(!"fun"%in%names(cols.args)){
      stop("No function provided to establish the color palette, please provide it")
    }
    if(!"cols"%in%names(cols.args)){
      cols.args<-c(cols.args,list("cols"=c("blue","yellow","red")))
      warning("No colors provided to establish the color palette, default blue-yellow-red color scale applied")
    }
    if(cols.args$fun=="discrete.palette"){
      DP.args<-cols.args[-which("fun"%in%names(cols.args))]
      if(!"ncols"%in%names(DP.args)){
        if(!is.null(mapping.args)&&"res"%in%names(mapping.args)){
          DP.args<-c(DP.args,list("ncols"=mapping.args$res+1))
        }
        else if(branch){
          DP.args<-c(DP.args,list("ncols"=length(DP.args$cols)))
        }
        else{
          DP.args<-c(DP.args,list("ncols"=1001))
        }
      }
      cols<-do.call("discrete.palette",DP.args)
    }
    else if(cols.args$fun=="scale.palette"){
      SP.args<-cols.args[-which("fun"%in%names(cols.args))]
      if(!"ncols"%in%names(SP.args)){
        if(!is.null(mapping.args)&&"res"%in%names(mapping.args)){
          SP.args<-c(SP.args,list("ncols"=mapping.args$res+1))
        }
        else if(branch){
          SP.args<-c(SP.args,list("ncols"=length(SP.args$cols)))
        }
        else{
          SP.args<-c(SP.args,list("ncols"=1001))
        }
      }
      if(!"middle.col"%in%names(SP.args)){SP.args<-c(SP.args,list("middle.col"=NA))}
      if(!"middle"%in%names(SP.args)){SP.args<-c(SP.args,list("middle"=NA))}
      if(!"span"%in%names(SP.args)){SP.args<-c(SP.args,list("span"=lims))}
      cols<-do.call("scale.palette",SP.args)
    }
  }
  else{
    cols<-cols.args
  }

  if(branch){
    col.scale<-TRUE
    if(!(length(order)==1&&order=="edge")){
      edge.order<-tree$edge[,2]
      edge.order[edge.order>Ntip(tree)]<-edge.order[edge.order>Ntip(tree)]-1
    }
    if(length(order)==1){
      if(order=="edge"&&length(cols)==length(values)){
        col.scale<-FALSE
      }
      if(order=="phylo"){
        values<-values[edge.order]
      }
      else if(order=="names"){
        values<-values[match(c(tree$tip.label,if(!is.null(tree$node.label)){tree$node.label[-1]}else{as.character((2:Nnode(tree))+Ntip(tree))}),names(values))][edge.order]
      }
    }
    else if(length(order)==length(values)){
      if(is.numeric(order)){
        values<-values[order(order)]
      }
      else if(is.character(order)){
        values<-values[match(c(tree$tip.label,if(!is.null(tree$node.label)){tree$node.label[-1]}else{as.character((2:Nnode(tree))+Ntip(tree))}),order)][edge.order]
      }
    }

    if(col.scale){
      cols<-cols[round(scales::rescale(values,c(1,length(cols)),lims))]
    }
    PP.args<-c(PP.args,list(x=tree,edge.color=cols))
    if(!"show.tip.label"%in%names(PP.args)){PP.args<-c(PP.args,list(show.tip.label=FALSE))}
    if(!"edge.width"%in%names(PP.args)){PP.args<-c(PP.args,list(edge.width=3))}
    par(mar=c(0,0,if(is.null(title)){0}else{2.1},0))
    do.call("plot.phylo",PP.args)
  }
  else{
    if(length(order)==length(values)){
      if(is.numeric(order)){
        values<-values[order(order)]
      }
      else if(is.character(order)){
        values<-values[if(!is.null(tree$node.label)){match(c(tree$tip.label,tree$node.label),order)}else{c(match(tree$tip.label,order[order%in%tree$tip.label]),which(!order%in%tree$tip.label))}]
      }
      order<-"phylo"
    }
    if(length(order)==1){
      if(order=="names"){
        tips<-values[tree$tip.label]
        if(!is.null(tree$node.label)){
          ancs<-unname(values[tree$node.label])
        }
        else{
          ancs<-unname(values[!names(values)%in%tree$tip.label])
        }
      }
      else if(order=="phylo"){
        tips<-setNames(values[1:Ntip(tree)],tree$tip.label)
        ancs<-unname(values[(1:Nnode(tree))+Ntip(tree)])
      }
      else{
        if(length(values)<(Ntip(tree)+Nnode(tree))){
          stop("Please provide as many values as there are taxa for the mapping")
        }
        else{
          tips<-setNames(values[-1][(tree$edge[,2])<=Ntip(tree)],tree$tip.label)
          ancs<-unname(c(values[1],values[-1][(tree$edge[,2])>Ntip(tree)]))
        }
      }
    }

    CM.args<-c(CM.args,list(tree=tree,x=tips,method="user",anc.states=ancs,plot=FALSE,lims=lims))
    CM<-do.call("contMap",CM.args)
    CM$cols<-replace(CM$cols,names(CM$cols),cols)

    PCM.args<-c(PCM.args,list(x=CM))
    if(!"legend"%in%names(PCM.args)){PCM.args<-c(PCM.args,list(legend=FALSE))}
    mar<-c(0,0,if(is.null(title)){0}else{2.1},0)
    if(!"mar"%in%names(PCM.args)){PCM.args<-c(PCM.args,list(mar=mar))}
    do.call("plot.contMap",PCM.args)
  }
  if(!is.null(title)){
    title(title)
  }
}
