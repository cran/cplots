#' @title Multi-class Stacked Circular Histogram and Rose Diagram
#' 
#' @description Function \code{cmhist} can be used to plot 2-dimensional 
#'   histograms and rose diagrams for circular data with multiple classes. 
#'   The histograms are stacked to avoid any overlap. 
#' 
#' @param value a numeric vector storing angular values between 0 and 2 pi,
#'   or an object that can be coerced to.
#' @param class a character vector specifying the group the \code{value}
#'   belongs to. It needs to have the same length as \code{value},
#'   otherwise it is repeated to the length of \code{value}. The order of
#'   plotting from the innermost to the outermost depends on the order of
#'   their appearance in \code{class}.
#' @param nbins the number of bins of the circular histogram. Internally,
#'   it is rounded to a multiple of 4.
#' @param radius the radius of the reference circle.  If \code{radius = 0},
#'   a rose diagram is produced; if \code{radius > 0}, a circular histogram
#'   is produced outside the reference circle.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param prob logical; if \code{TRUE}, the circular histogram graphic is a
#'   representation of probability densities; if \code{FALSE}, a
#'   representation of frequencies.
#' @param proportion logical; if \code{TRUE}, the frequencies are scaled by
#'   the proportion of each class, so that the total area under bars is
#'   unity; if \code{FALSE}, each class is considered as a separate
#'   distribution and has area of unity.
#' @param total.area a positive number specifying the total area under all the
#'   histograms. If \code{total.area = NULL}, no scaling is applied, the
#'   plot is in the original scale. If \code{area.prop = TRUE}, the total area 
#'   is automatically unity without scaling.
#' @param nlabels integer, for the number of levels to be plotted; if
#'   \code{0}, no label is plotted.  The larger the number is, the more
#'   accurate the plot will be.
#' @param cols the colors to fill the bars, with the same order as the
#'   class.
#' @param borders the colors of the border around the bars.
#' @param m the number of points within each bin to plot the circular
#'   histogram. The larger the number is, the smoother the plot looks.
#' @param xlim numeric vectors of length 2, giving the x coordinates
#'   ranges.
#' @param ylim numeric vectors of length 2, giving the y coordinates
#'   ranges.
#' @param main the main title (on top)
#' @param type the type of circular data, one of the values \code{"null"},
#'   \code{"compass"} or \code{"clock"}.  If \code{"null"}, no special
#'   lables plotted for directions. If \code{"compass"}, the four cardinal
#'   directions are printed inside the reference circle. If \code{"clock"},
#'   labels for 24 hours are printed inside the reference circle.
#' @param x.legend x coordinate to plot the legend.
#' @param y.legend y coordinate to plot the legend.
#' 
#' @concept multi-class stacked circular histogram
#' 
#' @return No return value
#' 
#' @author Danli Xu <dxu452@aucklanduni.ac.nz>, Yong Wang <yongwang@auckland.ac.nz>
#' 
#' @references
#'
#' Xu, D. and Wang, Y. (2020). Area-proportional Visualization for
#' Circular Data. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{29}, 351-357.
#' 
#' @seealso \code{\link{chist}}, \code{\link{cmdensity}}
#' 
#' @importFrom graphics hist plot points text legend polygon segments title
#' @importFrom stats uniroot
#' @importFrom grDevices hcl
#' 
#' @export
#' 
#' @examples
#' # Load the dataset
#' library(circular)
#' data("pigeons", package = "circular")
#' x = pigeons[,2] / 180 * pi
#' y = pigeons[,1]
#' 
#' # stacked circular histograms
#' cmhist(x, y)             # area-proportional 
#' cmhist(x, y, area=FALSE) # height-proportional
#' 


cmhist = function(value, class, nbins=36, radius=1/sqrt(base::pi),
                  area.prop=TRUE, prob=TRUE, proportion=FALSE, total.area=1, 
                  nlabels=4, cols=NULL, borders=NULL, m=NA, 
                  xlim=NULL, ylim=NULL, main=NULL,
                  type=c("null","compass","clock"),
                  x.legend="bottomright", y.legend=NULL) {
  value = as.vector(value)
  type = match.arg(type)
  pi = base::pi
  n = length(value)
  nbins = max(4, round(nbins / 4) * 4)
  nlabels = max(0, ceiling(nlabels))
  if(is.na(m)) m = max(ceiling(360 / nbins), 2)
  circle = seq(0, 2 * pi, len = 500)
  br = seq(0, 2 * pi, len = nbins + 1)
  br2 = br[-(nbins + 1)]  # remove the last element
  cb2 = cos(br2)
  sb2 = sin(br2)
  
  ## deal with multiple classes
  class = rep(class, length.out = n)
  level = factor(class, unique(class))              # preserve the order of class
  valuelist = split(value, level)                   # list of classified value
  prop = sapply(valuelist, length) / n              # proportion of each class
  k = length(prop)                                  # number of class
  histct =
    function(x, breaks = br) hist(x, breaks = breaks, plot = FALSE)$density
  ctmat = sapply(valuelist, histct)                 # raw density
  if(proportion) ctmat = ctmat * rep(prop*k, rep(nbins, k))  # scaled by prop.
  ctcum = t(apply(ctmat, 1, cumsum))                # cumsum
  ctmax = ctcum[, k]
  # factor = if(scale) scalefactor(ctmax, radius, k) else 1
  if (is.null(total.area)) factor = 1 
  else factor = scalefactor(ctmax, radius, total.area, area.prop) / 
         (if(area.prop) k else 1)
  ctcumf = circtrans(ctcum, radius, area.prop, factor)
  ctmaxf = ctcumf[,k]
    
  ## every point in the same nbins has the same height
  m1 = matrix(1:(nbins * m), nbins, m, byrow = TRUE)
  m2 = seq(m + 1, nbins * m + 1, m)
  mat = cbind(m1, m2)
  ind = t(mat)
  dim(ind) = NULL
  angle = seq(0, 2 * pi, len = nbins * m + 1)
  new = angle[ind]
  
  ## plot
  if (is.null(cols)) cols = hcl(seq(0, 240, len = k), c = 90, l = 80)
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(cb2 * ctmaxf))
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(sb2 * ctmaxf))
  plot(0, type="n", asp=1, bty="n", axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
  points(0, 0, pch = 3)
  switch(type,
         compass = 
           text(cos(0:3*pi/2) * radius * 0.8, 
                sin(0:3*pi/2) * radius * 0.8, 
                labels=c("E", "N", "W", "S"),
                cex = 1.2),
         clock = {
           timeseq = seq(0, 2 * pi, len = 25)[-25]
           text(cos(timeseq) * radius * 0.85, sin(timeseq) * radius * 0.85, 
                labels=c(6:1, 24:7))
         },
         null =, ) 
  
  ## labels
  if(nlabels > 0) {
    factor.fr = if(prob) 1 else 2 * pi * n / nbins
    ma = max(ctmax) * 1.1 * factor.fr
    by = if(prob) signif(ma/nlabels, 1) else max(signif(ma/nlabels,1), 1)
    digit = floor(log10(by))
    label = seq(0, ma, by=by)[-1]
    flabel = label / factor.fr             # area: density value for label
    lab = circtrans(flabel, radius, area.prop, factor)
    for (i in 1:length(label)) 
      lines(cos(circle) * lab[i], sin(circle) * lab[i], lty = 2, col="darkgrey")
    text(-lab/sqrt(2), lab/sqrt(2), cex=1.2, 
         labels=sprintf(paste0("%.", pmax(-digit, 0), "f"), label))
  }

  ## polygon set up
  for (i in 1:k) {
    xpoly = cos(new) * rep(ctcumf[,k-i+1], each = m + 1)
    ypoly = sin(new) * rep(ctcumf[,k-i+1], each = m + 1)
    circle2 = rev(circle)
    polygon(c(xpoly, cos(circle2) * radius), c(ypoly, sin(circle2) * radius), 
            col = cols[k:1][i], border = borders[k:1][i])
  }
  segments(x0=cb2*radius, y0=sb2*radius, x1=cb2*ctmaxf, y1=sb2*ctmaxf)
  
  ## main title
  if(is.null(main)) {
    main = paste0(if(area.prop) "Area" else "Height", "-proportional ",
                  if(radius != 0) "Circular Histogram" else "Rose Diagram")
  }
  title(main = main)

  if(!is.null(x.legend)) 
    legend(x.legend, y.legend, legend=levels(level), pch=22, pt.cex=2.2, pt.bg=cols,
           xjust=0.5, yjust=0.5, bg="white")
}





