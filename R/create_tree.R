#' @title Create a phylogenetic tree
#'
#' @description
#' This function creates a tree, more a less stepwise, depending on the data in the function.
#' The function can be run without arguments.
#' If the desired phylogeny is already known and prepared, it is advised to prepare and save the elements used in this function, in case of small errors (yet, the function doesn't allow one to recover inputted data if there is an error).
#'
#' @param nbtaxa Desired number of taxa; optional, especially if either \code{taxa} or \code{age_taxa} parameters are present
#' @param taxa A vector specifying the name of \strong{all} desired taxa; if not all taxa names are known, rather do not specify the known names yet and specify them by running the function without this parameter
#' @param age_taxa A vector specifying the age of \strong{all} desired taxa; if not all taxa ages are known, rather do not specify the known ages yet and specify them by running the function without this parameter
#' @param nbnodes Desired number of nodes; optional, especially if \code{age_nodes} is present
#' @param nodes A list of vectors; each vector is a node and must contains the names of the "targets", either the taxa (full name) or other nodes (Ni, i being the i-th node). Each list can be named as Ni (i being the i-th node) or not.
#' @param age_nodes A vector specifying the age of \strong{all} nodes; if not all node ages are known, rather do not specify the known ages yet and specify them by running the function without this parameter
#' @param format The format of the output; can be an object of class \code{phylo} by specifying \code{"phylo"} or \code{"phylo object"} or simply a newick/parenthetic text by specifying \code{"newick"}, \code{"NEWICK"} or \code{"parenthetic"}
#'
#' @return Returns either an object of class \code{phylo} or a text in newick/parenthetic format. If the latter, the tree can also be saved during performing the function.
#'
#' @importFrom graphics points text
#' @importFrom gtools ask
#' @importFrom ape plot.phylo read.tree write.tree
#' @importFrom foreach foreach %do%
#' @importFrom DescTools RoundTo
#'
#' @export

create.tree <- function(nbtaxa,taxa,age_taxa,nbnodes,nodes,age_nodes,format) {
  if(missing(age_taxa)){age_taxa<-NULL}
  if(missing(age_nodes)){age_nodes<-NULL}
  if(is.null(age_taxa)&is.null(age_nodes)){
    aged_tree<-ask("Do you want to give ages to the tips of the tree ? (Y/N)")
    if(aged_tree=="Y"|aged_tree=="y"|aged_tree=="YES"|aged_tree=="yes"|aged_tree=="Yes"){
      aged_tree<-TRUE
    }
    else{
      aged_tree<-FALSE
    }
  }
  else{
    aged_tree<-TRUE
  }
  ultra<-FALSE
  if(aged_tree==FALSE){
    ultra<-TRUE
  }
  else{
    if(is.null(age_nodes)){
      ultra<-ask("Do you want an ultrametric tree (i.e. without specifying ages of nodes) ? (Y/N)")
      if(ultra=="yes"|ultra=="Y"|ultra=="y"|ultra=="Yes"|ultra=="YES"){
        ultra<-TRUE
      }
    }
  }

  if (missing(taxa)) {
    taxa <- c()
    if (missing(nbtaxa)) {
      loop <- "go"
      nbtaxa <- 0
      while (loop != "end") {
        nbtaxa <- nbtaxa + 1
        temp_taxa <- ask(paste0("Please give the name of the ",nbtaxa,
                                if (nbtaxa == 1) {"-st"}
                                else{
                                  if (nbtaxa == 2) {"-nd"}
                                  else{"-th"}},
                                " taxon; if you finished, please write 'end'"))
        if (temp_taxa != "end") {
          temp_age_taxa <- ask(paste0("Please give the numeric age of the ",nbtaxa,
                                     if (nbtaxa == 1) {"-st"}
                                     else{
                                       if (nbtaxa == 2) {"-nd"}
                                       else{"-th"}},
                                     " taxon; if you finished, please write 'end'"))
          if(aged_tree){
            if (is.na(suppressWarnings(as.numeric(temp_age_taxa))) == TRUE) {
              cat("You did not provide an age for this taxon, please provide it now or function will stop","\n")
              temp_age_taxa <- ask(paste0("Please give the numeric age of the ",nbtaxa,
                                          if (nbtaxa == 1) {"-st"}
                                          else{
                                            if (nbtaxa == 2) {"-nd"}
                                            else{"-th"}},
                                          " taxon; if you finished, please write 'end'"))
            }
          }
          taxa[nbtaxa] <- as.character(temp_taxa)
          if(aged_tree){age_taxa[nbtaxa] <- as.numeric(temp_age_taxa)}
        }
        else{
          loop <- "end"
          nbtaxa <- nbtaxa - 1
        }
      }
    }
    else{
      for (i in 1:nbtaxa) {
        temp_taxa <- ask(paste0("Please give the name of the ", i,
                                if (i == 1) {"-st"}
                                else{
                                  if (i == 2) {"-nd"}
                                  else{"-th"}},
                                " taxon"))
        temp_age_taxa <- ask(paste0("Please give the numeric age of the ",i,
                                    if (i == 1) {"-st"}
                                    else{
                                      if (i == 2) {"-nd"}
                                      else{"-th"}},
                                    " taxon"))
        if(aged_tree){
          if (is.na(suppressWarnings(as.numeric(temp_age_taxa))) == TRUE) {
            cat("You did not provide an age for this taxon, please provide it now or function will stop","\n")
            temp_age_taxa <- ask(paste0("Please give the numeric age of the ",i,
                                        if (i == 1) {"-st"}
                                        else{
                                          if (i == 2) {"-nd"}
                                          else{"-th"}},
                                        " taxon; if you finished, please write 'end'"))
          }
        }
        taxa[i] <- as.character(temp_taxa)
        if(aged_tree){age_taxa[i] <- as.numeric(temp_age_taxa)}
      }
    }
  }
  else{
    taxa <- as.character(taxa)
    nbtaxa <- length(taxa)
    if(aged_tree){
      if (all(suppressWarnings(is.na(as.numeric(age_taxa))))) {
        for (i in 1:length(taxa)) {
          temp_age_taxa <- ask(paste0("Please give the numeric age of the taxon ", taxa[i]))
          if (all(is.na(suppressWarnings(as.numeric(temp_age_taxa))) == TRUE)) {
            cat("You did not provide an age for this taxon, please provide it now or function will stop","\n")
            temp_age_taxa <- ask(paste0("Please give the numeric age of the taxon ", taxa[i]))
          }
          age_taxa[i] <- as.numeric(temp_age_taxa)
        }
        age_taxa <- as.numeric(age_taxa)
      }
      else{
        age_taxa <- as.numeric(age_taxa)
      }
    }
  }
  if(aged_tree==FALSE){
    age_taxa<-rep(0,length(taxa))
  }
  if(missing(nodes)==TRUE){
    temp_phylo <- c()
    plot_depth <- max(age_taxa,if(is.null(age_nodes)==FALSE){age_nodes})
    for (i in 1:nbtaxa) {
      temp_phylo[i] <- paste0(taxa[i], ":", plot_depth - age_taxa[i], ",", collapse = "")
    }
    temp_phylo <- paste0("(", gsub('.{1}$', '', paste(temp_phylo, collapse = "", sep = "")), ");")
    temp <- read.tree(text = temp_phylo)
    taxa_short<-c()
    for (i in 1:length(taxa)){
      if(nchar(taxa[i])>4){
        if(length(strsplit(taxa[i],split="[ ]")[[1]])>1){
          temp_tax_names<-strsplit(taxa[i],split="[ ]")[[1]]
          taxa_short[i]<-paste0(unlist(foreach(j=1:length(temp_tax_names))%do%paste0(strsplit(temp_tax_names,split="")[[j]][1:2],collapse="")),collapse=" ")
        }
        else{
          taxa_short[i]<-paste0(strsplit(taxa[i],split="")[[1]][1:4],collapse="")
        }
      }
      else{
        taxa_short[i]<-taxa[i]
      }
    }
    temp$tip.label <- taxa_short
    x_gap<-(plot_depth - min(age_taxa))/10
    if(x_gap==0){x_gap<-1}
    if(log10(x_gap)+1<1){
      x_gap_magnitude<-0
      x_gap_temp<-0
      while(x_gap_temp<1){
        x_gap_temp<-x_gap*10^x_gap_magnitude
        if(x_gap_temp<1){x_gap_magnitude<-x_gap_magnitude+1}
      }
      x_gap<-RoundTo(round(x_gap,x_gap_magnitude),5/10^x_gap_magnitude)
    }
    else{
      x_gap<-ceiling(x_gap)
    }
    plot.phylo(temp,edge.color = "white",x.lim = c(min(age_taxa), plot_depth),label.offset = -x_gap/10)
    x_ticks <-invisible(axis(1,labels = NA,tick = FALSE,at = seq(floor(min(age_taxa)), plot_depth, x_gap)+(RoundTo(plot_depth,x_gap)-plot_depth)))
    x_labels <- rev(x_ticks)-(RoundTo(plot_depth,x_gap)-plot_depth)
    axis(1, at = x_ticks, labels = x_labels)
    node_nb_ends <- list()
    node_ends <- list()
    node_coords_x <- c()
    node_coords_y <- c()
    loop <- "go"
    node <- 0
    nodemax<--1
    if(is.null(age_nodes)==FALSE|missing(nbnodes)==FALSE){
      if(missing(nbnodes)){
        nbnodes<-length(age_nodes)
      }
      nodemax<-nbnodes
    }
    else{
      nbnodes<-c()
    }
    while (loop != "end" & node != nodemax) {
      node <- node + 1
      node_nb_ends[[node]] <- ask(paste0(
        if(nodemax<0){"If their are no more nodes, type 'end', otherwise, i"}
        else{"I"},
        "f the ",node,
        if (node == 1) {"-st"}
        else{
          if (node == 2) {"-nd"}
          else{"-th"}},
        " node is dichotomic, type 2, or type the number of branches arising from it"))
      if (node_nb_ends[[node]] == "end") {
        loop <- "end"
        node_nb_ends <- node_nb_ends[-node]
        node <- node - 1
        nbnodes<-node
        break
      }
      if (is.na(suppressWarnings(as.numeric(node_nb_ends[[node]]))) == TRUE) {
        cat("You did not provide a number of branches, please provide it now or function will stop","\n")
        node_nb_ends[[node]] <- ask(paste0(
          if(nodemax<0){"If their are no more nodes, type 'end', otherwise, i"}
          else{"I"},
          "f the ",node,
          if (node == 1) {"-st"}
          else{
            if (node == 2) {"-nd"}
            else{"-th"}},
          " node is dichotomic, type 2, or type the number of branches arising from it"))
      }
      node_ends[[node]] <- as.character(0)
      for (i in 1:as.numeric(node_nb_ends[[node]])) {
        node_ends[[node]][i] <- ask(paste0("If the ",i,
                                           if (i == 1) {"-st"}
                                           else{
                                             if (i == 2) {"-nd"}
                                             else{"-th"}},
                                           " end is a taxon, type its name; if it is a node, type 'N' and its displayed number"))
        if (((any(taxa == node_ends[[node]][i])) | (any(names(node_ends) == node_ends[[node]][i]))) == FALSE) {
          cat("You did not provide the name of the written taxa or nodes, please provide it now or function will stop","\n")
          node_ends[[node]][i] <- ask(paste0("If the ",i,
                                             if (i == 1) {"-st"}
                                             else{
                                               if (i == 2) {"-nd"}
                                               else{"-th"}},
                                             " end is a taxon, type its name; if it is a node, type 'N' and its displayed number"))
        }
        if (suppressWarnings(is.na(as.numeric(paste(strsplit(node_ends[[node]][i], "")[[1]][-1],collapse=""))))==TRUE) {
          node_ends[[node]][i] <- which(taxa == node_ends[[node]][i])
        }
        names(node_ends)[node] <- paste0("N", node)
      }
      if(nodemax>0){
        age_nodes<-age_nodes
      }
      else{
        if(ultra){
          age_nodes[node] <- max(c(suppressWarnings(max(age_taxa[na.omit(as.numeric(node_ends[[node]]))])), suppressWarnings(max(age_nodes[unlist(foreach(i=1:as.numeric(node_nb_ends[[node]]))%do%which((names(node_ends)==node_ends[[node]][i])))])))) + 1
        }
        else{
          age_nodes[node] <- as.numeric(ask(paste0("Please give the age of the ",node,
                                                   if (node == 1) {"-st"}
                                                   else{
                                                     if (node == 2) {"-nd"}
                                                     else{"-th"}},
                                                   " node; if not known, type 0")))
        }
        if (age_nodes[node] < suppressWarnings(max(age_taxa[na.omit(as.numeric(node_ends[[node]]))]))) {
          age_nodes[node] <- as.numeric(ask(paste0("Error: node age is more recent than (at least) on of the terminal taxa; please give its age or, if not known, type 0")))
        }
        if (age_nodes[node] == 0) {
          if (length(age_nodes) == 1) {
            age_nodes[node] <-suppressWarnings(max(age_taxa[na.omit(as.numeric(node_ends[[node]]))])) + 1
          }
          else{
            age_nodes[node] <- max(c(suppressWarnings(max(age_taxa[na.omit(as.numeric(node_ends[[node]]))])), suppressWarnings(max(age_nodes[unlist(foreach(i=1:as.numeric(node_nb_ends[[node]]))%do%which((names(node_ends)==node_ends[[node]][i])))])))) + 1
          }
        }
      }
      max_depth <- max(c(max(age_taxa), max(age_nodes)),na.rm=TRUE)
      node_coords_x[node] <- plot_depth - age_nodes[node]
      if (all(suppressWarnings(is.na(as.numeric(node_ends[[node]]))) == FALSE)) {
        node_coords_y[node] <- mean(as.numeric(node_ends[[node]]))
      }
      else{
        where <- suppressWarnings(is.na(as.numeric(node_ends[[node]])))
        if (length(levels(as.factor(where))) == 2) {
          prev_node <- as.numeric(paste(strsplit(node_ends[[node]][where == TRUE], "")[[1]][-1],collapse=""))
          node_coords_y[node] <- mean(c(as.numeric(node_ends[[node]][where == FALSE]), node_coords_y[prev_node]))
        }
        else{
          prev_node_1 <- as.numeric(paste(strsplit(node_ends[[node]][where == TRUE], "")[[1]][-1],collapse=""))
          prev_node_2 <- as.numeric(paste(strsplit(node_ends[[node]][where == TRUE], "")[[2]][-1],collapse=""))
          node_coords_y[node] <- mean(c(node_coords_y[prev_node_1], node_coords_y[prev_node_2]))
        }
      }
      if (max_depth > plot_depth) {
        plot.phylo(temp,edge.color = "white",x.lim = c((min(age_taxa) - (max_depth - max(age_taxa))), max(age_taxa)),label.offset = 0)
        x_ticks <- invisible(axis(1,labels = NA,tick = FALSE,at = seq(floor(min(age_taxa) - (max_depth - max(age_taxa))), ceiling(max(age_taxa)), max(1,round(log10(max(age_taxa) - (min(age_taxa) - (max_depth - max(age_taxa)))))))))
        if (min(x_ticks) < 0) {
          x_labels <- rev(x_ticks - (min(x_ticks)))
        }
        else{
          x_labels <- rev(x_ticks)
        }
        axis(1, at = x_ticks, labels = x_labels)
        for (i in 1:max(c(1, length(node_ends) - 1))) {
          text(node_coords_x[i],node_coords_y[i],labels = names(node_ends)[i],offset = 0)
          x_offset <- strwidth("i")
          for (j in 1:as.numeric(node_nb_ends[[i]])) {
            start_x <- node_coords_x[i] + x_offset
            end_x <- c()
            if (suppressWarnings(is.na(as.numeric(node_ends[[i]][j]))) == FALSE) {
              end_x <- node_coords_x[i] + age_nodes[i] - age_taxa[as.numeric(node_ends[[i]][j])] - x_offset
            }
            else{
              end_x <- node_coords_x[i] + age_nodes[i] - age_nodes[as.numeric(paste(strsplit(node_ends[[i]][j], split = "")[[1]][-1],collapse=""))] - x_offset
            }
            start_y <- node_coords_y[i]
            end_y <- c()
            if (suppressWarnings(is.na(as.numeric(node_ends[[i]][j]))) == FALSE) {
              end_y <- as.numeric(node_ends[[i]][j])
            }
            else{
              end_y <- node_coords_y[as.numeric(paste(strsplit(node_ends[[i]][j], split = "")[[1]][-1],collapse=""))]
            }
            points(c(start_x, end_x), c(start_y, end_y), type = "l")
          }
        }
      }
      text(node_coords_x[node],node_coords_y[node],labels = names(node_ends)[node],offset = 0)
      x_offset <- strwidth("i")
      for (i in 1:as.numeric(node_nb_ends[[node]])) {
        start_x <- node_coords_x[node] + x_offset
        end_x <- c()
        if (suppressWarnings(is.na(as.numeric(node_ends[[node]][i]))) == FALSE) {
          end_x <- node_coords_x[node] + age_nodes[node] - age_taxa[as.numeric(node_ends[[node]][i])] - x_offset
        }
        else{
          end_x <- node_coords_x[node] + age_nodes[node] - age_nodes[as.numeric(paste(strsplit(node_ends[[node]][i], split = "")[[1]][-1],collapse=""))] - x_offset
        }
        start_y <- node_coords_y[node]
        end_y <- c()
        if (suppressWarnings(is.na(as.numeric(node_ends[[node]][i]))) == FALSE) {
          end_y <- as.numeric(node_ends[[node]][i])
        }
        else{
          end_y <- node_coords_y[as.numeric(paste(strsplit(node_ends[[node]][i], split = "")[[1]][-1],collapse=""))]
        }
        points(c(start_x, end_x), c(start_y, end_y), type = "l")
      }
    }
  }
  else{
    nbnodes<-length(nodes)
    node_nb_ends<-foreach(i=1:length(nodes))%do%length(nodes[[i]])
    node_ends<-nodes
    if(is.null(names(node_ends))){
      names(node_ends)<-rep("NA",length(node_ends))
    }
    for (i in 1:nbnodes){
      for (j in 1:node_nb_ends[[i]]){
        if (suppressWarnings(is.na(as.numeric(paste(strsplit(node_ends[[i]][j], "")[[1]][-1],collapse=""))))==TRUE) {
          node_ends[[i]][j] <- which(taxa == node_ends[[i]][j])
        }
      }
      if(paste0("N",i)!=names(node_ends)[[i]]){
        names(node_ends)[[i]]<-paste0("N",i)
      }
      if(ultra){
        age_nodes[i] <- max(c(suppressWarnings(max(age_taxa[na.omit(as.numeric(node_ends[[i]]))])), suppressWarnings(max(age_nodes[unlist(foreach(j=1:as.numeric(node_nb_ends[[i]]))%do%which((names(node_ends)==node_ends[[i]][j])))])))) + 1
      }
    }
  }

  numb_branches <- lengths(node_ends)
  branches_temp <- matrix(ncol = 2 * max(numb_branches),nrow = length(node_ends),NA)
  for (i in 1:length(node_ends)) {
    term_unit <- c(node_ends[[i]])
    temp <- c()
    for (j in 1:numb_branches[i]) {
      if (is.na(suppressWarnings(as.numeric(term_unit[j]))) == TRUE) {
        temp[j] <- as.numeric(paste(strsplit(term_unit[j], "")[[1]][-1],collapse="")) + length(taxa)
      }
      else{
        temp[j] <- as.numeric(term_unit[j])
      }
    }
    term_unit <- temp
    apic_unit <- as.numeric(paste(strsplit(names(node_ends)[i], "")[[1]][-1],collapse="")) + length(taxa)
    branches_temp[i, seq(1, 1 + 2 * (numb_branches[i] - 1), 2)] <- apic_unit
    branches_temp[i, seq(2, 2 + 2 * (numb_branches[i] - 1), 2)] <- term_unit
  }
  branches <- na.omit(matrix(ncol = 2,nrow = prod(dim(branches_temp)) / 2,t(branches_temp),byrow = TRUE))
  for (i in 1:2) {
    branches[, i] <- rev(branches[, i])
  }
  age_all <- c(age_taxa, age_nodes)
  age_branches <- c()
  for (i in 1:length(branches[, 1])) {
    age_branches[i] <- age_all[branches[i, 1]] - age_all[branches[i, 2]]
  }
  if(aged_tree){
    output <-list(edge = branches,edge.length = age_branches,Nnode = as.integer(nbnodes),tip.label = taxa)
  }
  else{
    output <-list(edge = branches,Nnode = as.integer(nbnodes),tip.label = taxa)
  }
  attributes(output)$class <- "phylo"
  attributes(output)$order <- "cladewise"
  output <- read.tree(text = write.tree(output))
  if(suppressWarnings(is.null(plot.phylo(output)))==FALSE){
    if (missing(format)) {
      format <- "phylo object"
    }
    if (format == "parenthetic" | format == "newick" | format == "NEWICK") {
      q_save <- ask("Do you want to save the tree ? (Y/N)")
      if (q_save == "Y" | q_save == "yes" | q_save == "Yes" | q_save == "YES" | q_save == "y") {
        loc <- paste0(
          ask("Please give the path where you want to save your tree"),
          ask("Please give the name of your tree"),
          ".tre"
        )
        write.tree(output, file = loc)
      }
      else{
        output <- write.tree(output)
      }
    }
    if (format == "phylo object" | format == "phylo"){
      output <- output
    }
    if(format=="")
    return(output)
  }
  else{
    output_error<-list(taxa=taxa,age_taxa=age_taxa,age_nodes=age_nodes,node_ends=node_ends,node_nb_ends=node_nb_ends,numb_branches=numb_branches,branches_temp=branches_temp,branches=branches,age_branches=age_branches)
    return(output_error)
  }
}
