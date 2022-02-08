#' @title Custom phytools:::plotPhylogram function to allow vertical space between taxa of a horizontal phylogram
#' @param x Mapping plot of phytools function (e.g., \code{plot.contMap}, \code{plotSimmap}) containing a phylogenetic \code{tree} and a vector of colors named \code{cols}
#' @param direction Horizontal direction of the phylogram. Either \code{"rightwards"} or \code{"leftwards"}
#' @param x_gap Numeric horizontal free space at the side of the plot, useful if one wants to contrast two trees. Default is 0
#' @param ylim Classic \code{ylim} vector. Default is sum of number of tips and of vertical spaces between taxa
#' @param y.label.offset Vector of same length than the number of tips of the tree. Needs to specify the offset of each tip (0 if same place, positive values otherwise).
#' @param fsize Relative font size for tip labels.
#' @param ftype Font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics)
#' @param lwd Line width for plotting
#' @param pts Logical value indicating whether or not to plot filled circles at each vertex of the tree, as well as at transition points between mapped states. Default is \code{FALSE}
#' @param node.numbers Logical value indicating whether or not node numbers should be plotted. Default is \code{FALSE}
#' @param mar Vector containing the margins for the plot to be passed to \code{\link{par}}. If not specified, the default margins all are 0
#' @param add Logical value indicating whether or not to add the plotted tree to the current plot (\code{TRUE}) or create a new plot (\code{FALSE}, the default)
#' @param offset Offset for the tip labels in character widths. Default is 0
#' @param setEnv Logical value indicating whether or not to set the environment \code{.PlotPhyloEnv}. Setting this to \code{TRUE} (the default) will allow compatibility with ape labeling functions such as \code{\link{nodelabels}}
#' @param placement Node placement following Felsenstein (2004; pp. 574-576). Can be \code{"intermediate"} (the default), \code{"centered"}, \code{"weighted"}, or \code{"inner"}
#' @param tips Labeled vector containing the vertical position of tips. Normally this will be \code{1:N} for \code{N} tips in the tree.
#' @param split.vertical Split the color of the vertically plotted edges by the state of the daughter edges. Only applies if the edge state changes exactly at a node
#' @param lend Line end style. See \code{\link{par}}
#' @param asp Aspect ratio. See \code{\link{plot.window}}
#' @param plot Logical value indicating whether or not to actually plot the tree. (See equivalent argument in \code{\link{plot.phylo}}
#' @details See \code{\link[phytools]{plotSimmap}} help page for more details
#' @keywords internal

plotPhylogram.JM<-function (x,
                            direction="rightwards",
                            x_gap=0,
                            ylim,
                            y.label.offset,
                            fsize=0.75,
                            ftype="off",
                            lwd=4,
                            pts=FALSE,
                            node.numbers=FALSE,
                            mar=rep(0,4),
                            add=FALSE,
                            offset=0,
                            setEnv=TRUE,
                            placement="intermediate",
                            tips=NULL,
                            split.vertical=FALSE,
                            lend=2,
                            asp=NA,
                            plot=TRUE) {
  tree<-x$tree
  colors<-x$cols
  ftype<-which(c("off", "reg", "b", "i", "bi") == ftype) - 1

  xlim<-if(ftype==0){if(direction=="leftwards"){c(-x_gap,max(nodeHeights(tree)))}else{c(0,max(nodeHeights(tree))+x_gap)}}else{xlim<-NULL}

  if (split.vertical && !setEnv) {
    cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
    spit.vertical <- FALSE
  }
  offsetFudge <- 1.37
  cw <- reorderSimmap(tree)
  pw <- reorderSimmap(tree, "postorder")
  n <- Ntip(cw)
  m <- cw$Nnode
  Y <- matrix(NA, m + n, 1)
  if(missing(y.label.offset)){y.label.offset<-rep(0,Ntip(tree))}else{y.label.offset<-cumsum(y.label.offset)}
  if(missing(ylim)){ylim<-c(1,Ntip(tree)+max(y.label.offset))}
  if (is.null(tips))
    Y[cw$edge[cw$edge[, 2] <= n, 2]] <- 1:n+y.label.offset
  else Y[cw$edge[cw$edge[, 2] <= n, 2]] <- if (is.null(names(tips)))
    tips[sapply(1:Ntip(cw), function(x, y) which(y == x),
                y = cw$edge[cw$edge[, 2] <= n, 2])]
  else tips[gsub("_", " ", cw$tip.label)]
  nodes <- unique(pw$edge[, 1])
  for (i in 1:m) {
    if (placement == "intermediate") {
      desc <- pw$edge[which(pw$edge[, 1] == nodes[i]),
                      2]
      Y[nodes[i]] <- (min(Y[desc]) + max(Y[desc]))/2
    }
    else if (placement == "centered") {
      desc <- getDescendants(pw, nodes[i])
      desc <- desc[desc <= Ntip(pw)]
      Y[nodes[i]] <- (min(Y[desc]) + max(Y[desc]))/2
    }
    else if (placement == "weighted") {
      desc <- pw$edge[which(pw$edge[, 1] == nodes[i]),
                      2]
      n1 <- desc[which(Y[desc] == min(Y[desc]))]
      n2 <- desc[which(Y[desc] == max(Y[desc]))]
      v1 <- pw$edge.length[which(pw$edge[, 2] == n1)]
      v2 <- pw$edge.length[which(pw$edge[, 2] == n2)]
      Y[nodes[i]] <- ((1/v1) * Y[n1] + (1/v2) * Y[n2])/(1/v1 +
                                                          1/v2)
    }
    else if (placement == "inner") {
      desc <- getDescendants(pw, nodes[i])
      desc <- desc[desc <= Ntip(pw)]
      mm <- which(abs(Y[desc] - median(Y[1:Ntip(pw)])) ==
                    min(abs(Y[desc] - median(Y[1:Ntip(pw)]))))
      if (length(mm > 1))
        mm <- mm[which(Y[desc][mm] == min(Y[desc][mm]))]
      Y[nodes[i]] <- Y[desc][mm]
    }
  }
  H <- nodeHeights(cw)
  par(mar = mar)
  if (is.null(offset))
    offset <- 0.2 * lwd/3 + 0.2/3
  if (!add)
    plot.new()
  if (is.null(xlim)) {
    pp <- par("pin")[1]
    sw <- fsize * (max(strwidth(cw$tip.label, units = "inches"))) +
      offsetFudge * fsize * strwidth("W", units = "inches")
    alp <- optimize(function(a, H, sw, pp) (a * 1.04 * max(H) +
                                              sw - pp)^2, H = H, sw = sw, pp = pp, interval = c(0,
                                                                                                1e+06))$minimum
    xlim <- if (direction == "leftwards")
      c(min(H) - sw/alp, max(H))
    else c(min(H), max(H) + sw/alp)
  }
  if (is.null(ylim))
    ylim = range(Y)
  if (direction == "leftwards")
    H <- max(H) - H
  plot.window(xlim = xlim, ylim = ylim, asp = asp)
  if (plot) {
    if (!split.vertical) {
      for (i in 1:m) lines(H[which(cw$edge[, 1] == nodes[i]),
                             1], Y[cw$edge[which(cw$edge[, 1] == nodes[i]),
                                           2]], col = colors[names(cw$maps[[match(nodes[i],
                                                                                  cw$edge[, 1])]])[1]], lwd = lwd)
    }
    for (i in 1:nrow(cw$edge)) {
      x <- H[i, 1]
      for (j in 1:length(cw$maps[[i]])) {
        if (direction == "leftwards")
          lines(c(x, x - cw$maps[[i]][j]), c(Y[cw$edge[i,
                                                       2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]],
                lwd = lwd, lend = lend)
        else lines(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i,
                                                          2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]],
                   lwd = lwd, lend = lend)
        if (pts)
          points(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i,
                                                        2]], Y[cw$edge[i, 2]]), pch = 20, lwd = (lwd -
                                                                                                   1))
        x <- x + if (direction == "leftwards")
          -cw$maps[[i]][j]
        else cw$maps[[i]][j]
        j <- j + 1
      }
    }
    if (node.numbers) {
      symbols(if (direction == "leftwards")
        max(H)
        else 0, mean(Y[cw$edge[cw$edge[, 1] == (Ntip(cw) +
                                                  1), 2]]), rectangles = matrix(c(1.2 * fsize *
                                                                                    strwidth(as.character(Ntip(cw) + 1)), 1.4 *
                                                                                    fsize * strheight(as.character(Ntip(cw) + 1))),
                                                                                1, 2), inches = FALSE, bg = "white", add = TRUE)
      text(if (direction == "leftwards")
        max(H)
        else 0, mean(Y[cw$edge[cw$edge[, 1] == (Ntip(cw) +
                                                  1), 2]]), Ntip(cw) + 1, cex = fsize)
      for (i in 1:nrow(cw$edge)) {
        x <- H[i, 2]
        if (cw$edge[i, 2] > Ntip(cw)) {
          symbols(x, Y[cw$edge[i, 2]], rectangles = matrix(c(1.2 *
                                                               fsize * strwidth(as.character(cw$edge[i,
                                                                                                     2])), 1.4 * fsize * strheight(as.character(cw$edge[i,
                                                                                                                                                        2]))), 1, 2), inches = FALSE, bg = "white",
                  add = TRUE)
          text(x, Y[cw$edge[i, 2]], cw$edge[i, 2], cex = fsize)
        }
      }
    }
    if (direction == "leftwards")
      pos <- if (par()$usr[1] > par()$usr[2])
        4
    else 2
    if (direction == "rightwards")
      pos <- if (par()$usr[1] > par()$usr[2])
        2
    else 4
    for (i in 1:n) if (ftype)
      text(H[which(cw$edge[, 2] == i), 2], Y[i], cw$tip.label[i],
           pos = pos, offset = offset, cex = fsize, font = ftype)
  }
  if (setEnv) {
    PP <- list(type = "phylogram", use.edge.length = TRUE,
               node.pos = 1, show.tip.label = if (ftype) TRUE else FALSE,
               show.node.label = FALSE, font = ftype, cex = fsize,
               adj = 0, srt = 0, no.margin = FALSE, label.offset = offset,
               x.lim = xlim, y.lim = ylim, direction = direction,
               tip.color = "black", Ntip = Ntip(cw), Nnode = cw$Nnode,
               edge = cw$edge, xx = sapply(1:(Ntip(cw) + cw$Nnode),
                                           function(x, y, z) y[match(x, z)], y = H, z = cw$edge),
               yy = Y[, 1])
    assign("last_plot.phylo", PP, envir = .PlotPhyloEnv)
  }
  if (plot)
    if (split.vertical)
      splitEdgeColor(cw, colors, lwd)
}
