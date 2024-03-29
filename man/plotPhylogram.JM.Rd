% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPhylogram_JM.R
\name{plotPhylogram.JM}
\alias{plotPhylogram.JM}
\title{Custom phytools:::plotPhylogram function to allow vertical space between taxa of a horizontal phylogram}
\usage{
plotPhylogram.JM(
  x,
  direction = "rightwards",
  x_gap = 0,
  ylim,
  y.label.offset,
  fsize = 0.75,
  ftype = "off",
  lwd = 4,
  pts = FALSE,
  node.numbers = FALSE,
  mar = rep(0, 4),
  add = FALSE,
  offset = 0,
  setEnv = TRUE,
  placement = "intermediate",
  tips = NULL,
  split.vertical = FALSE,
  lend = 2,
  asp = NA,
  plot = TRUE
)
}
\arguments{
\item{x}{Mapping plot of phytools function (e.g., \code{plot.contMap}, \code{plotSimmap}) containing a phylogenetic \code{tree} and a vector of colors named \code{cols}}

\item{direction}{Horizontal direction of the phylogram. Either \code{"rightwards"} or \code{"leftwards"}}

\item{x_gap}{Numeric horizontal free space at the side of the plot, useful if one wants to contrast two trees. Default is 0}

\item{ylim}{Classic \code{ylim} vector. Default is sum of number of tips and of vertical spaces between taxa}

\item{y.label.offset}{Vector of same length than the number of tips of the tree. Needs to specify the offset of each tip (0 if same place, positive values otherwise).}

\item{fsize}{Relative font size for tip labels.}

\item{ftype}{Font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics)}

\item{lwd}{Line width for plotting}

\item{pts}{Logical value indicating whether or not to plot filled circles at each vertex of the tree, as well as at transition points between mapped states. Default is \code{FALSE}}

\item{node.numbers}{Logical value indicating whether or not node numbers should be plotted. Default is \code{FALSE}}

\item{mar}{Vector containing the margins for the plot to be passed to \code{\link{par}}. If not specified, the default margins all are 0}

\item{add}{Logical value indicating whether or not to add the plotted tree to the current plot (\code{TRUE}) or create a new plot (\code{FALSE}, the default)}

\item{offset}{Offset for the tip labels in character widths. Default is 0}

\item{setEnv}{Logical value indicating whether or not to set the environment \code{.PlotPhyloEnv}. Setting this to \code{TRUE} (the default) will allow compatibility with ape labeling functions such as \code{\link{nodelabels}}}

\item{placement}{Node placement following Felsenstein (2004; pp. 574-576). Can be \code{"intermediate"} (the default), \code{"centered"}, \code{"weighted"}, or \code{"inner"}}

\item{tips}{Labeled vector containing the vertical position of tips. Normally this will be \code{1:N} for \code{N} tips in the tree.}

\item{split.vertical}{Split the color of the vertically plotted edges by the state of the daughter edges. Only applies if the edge state changes exactly at a node}

\item{lend}{Line end style. See \code{\link{par}}}

\item{asp}{Aspect ratio. See \code{\link{plot.window}}}

\item{plot}{Logical value indicating whether or not to actually plot the tree. (See equivalent argument in \code{\link{plot.phylo}}}
}
\description{
Custom phytools:::plotPhylogram function to allow vertical space between taxa of a horizontal phylogram
}
\details{
See \code{\link[phytools]{plotSimmap}} help page for more details
}
\keyword{internal}
