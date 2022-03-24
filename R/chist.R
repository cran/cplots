#' @title Circular Histogram and Rose Diagram
#' 
#' @description Function \code{chist} can be used to plot 2-dimensional 
#'   histograms and rose diagrams for circular data.
#' 
#' @param x a numeric vector storing angular values between 0 and 2 pi, or
#'   an object that can be coerced to.
#' @param nbins the number of bins of the circular histogram.  Internally,
#'   it is rounded to a multiple of 4.
#' @param radius the radius of the reference circle. If \code{radius = 0},
#'   a rose diagram is produced; if \code{radius > 0}, a circular histogram
#'   is produced outside the reference circle.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param prob logical; if \code{TRUE}, the circular histogram graphic is a
#'   representation of probability densities; if \code{FALSE}, a
#'   representation of frequencies.
#' @param total.area a positive number specifying the total area under the
#'   density curve. If \code{total.area = NULL}, no scaling is applied, the
#'   plot is in the original scale. If \code{area.prop = TRUE}, the total area 
#'   is automatically unity without scaling.
#' @param nlabels integer, for the number of levels for the
#'   density/frequency values to be plotted; if \code{0}, no label is
#'   plotted
#' @param col the color to fill the bars.
#' @param border the color of the border around the bars.
#' @param m the number of points within each bin to plot the circular
#'   histogram. The larger the number is, the smoother the plot looks.
#' @param xlim numeric vectors of length 2, giving the x coordinates
#'   ranges.
#' @param ylim numeric vectors of length 2, giving the y coordinates
#'   ranges.
#' @param main the main title (on top)
#' 
#' @concept circular histogram
#' @concept rose diagram
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
#' @seealso \code{\link{cbarplot}}, \code{\link{cdensity}}, \code{\link{cdotplot}}
#' 
#' @importFrom graphics hist plot abline text segments points title
#' @importFrom stats uniroot
#' @importFrom circular rvonmises circular
#' 
#' @export
#' 
#' @examples
#' # 600 observations from two von Mises distributions
#' library(circular)
#' x = c(rvonmises(200, circular(pi/4), 5), rvonmises(400, circular(pi), 20))
#'
#' chist(x)                           # area-proportional circular histgram
#' chist(x, area = FALSE)             # height-proportional circular histgram
#' chist(x, radius=0)                 # area-proportional rose diagram
#' chist(x, radius=0, area=FALSE)     # height-proportional rose diagram
#' 
#' chist(x, prob=FALSE)               # labels for frequency
#' chist(x, nlabels=0)                # no label
#' chist(x, xlim=c(-1.7,1))           # use xlim
#' chist(x, area=FALSE, total=2)      # with scaling
#' chist(x, area=FALSE, total=NULL)   # without scaling
#' 


chist = function(x, nbins=36, radius=1/sqrt(base::pi), area.prop=TRUE,
                 prob=TRUE, total.area=1, nlabels=4,
                 col="lightblue", border="skyblue4", m=NA,  
                 xlim=NULL, ylim=NULL, main=NULL) {
  x = as.vector(x)
  n = length(x)
  pi = base::pi
  nbins = max(4, round(nbins / 4) * 4)
  nlabels = max(0, ceiling(nlabels))
  if(is.na(m)) m = max(ceiling(360 / nbins), 2)
  circle = seq(0, 2 * pi, len = 200)
  br = seq(0, 2 * pi, len = nbins + 1)                # break points
  d = hist(x, plot = FALSE, breaks = br)$density      # density
  # factor = if(scale) scalefactor(d, radius) else 1
  factor = scalefactor(d, radius, total.area, area.prop)
  if (is.null(total.area)) factor = 1 
  df = circtrans(d, radius, area.prop, factor)

  ## every point in the same bin has the same height
  m1 = matrix(1:(nbins * m), nbins, m, byrow = TRUE)
  m2 = seq(m + 1, nbins * m + 1, by = m)
  mat = cbind(m1, m2)
  ind = t(mat)
  dim(ind) = NULL
  angle = seq(0, 2 * pi, len = nbins * m + 1)
  new = angle[ind]
  cb = cos(br[-(nbins+1)])
  sb = sin(br[-(nbins+1)])
  x1 = cb * df
  y1 = sb * df
  
  ## plot
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(x1))
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(y1))
  plot(0, type="n", asp=1, bty="n", axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
  points(0, 0, pch = 3)
  if(nlabels > 0) {
    factor.d = if(prob) 1 else 2 * pi * n / nbins
    ma = max(d) * 1.1 * factor.d
    by = if(prob) signif(ma/nlabels, 1) else max(signif(ma/nlabels,1), 1)
    digit = floor(log10(by))
    label = seq(0, ma+by/2, by=by)[-1]
    flabel = label / factor.d # density values for labelling
    lab = circtrans(flabel, radius, area.prop, factor)
    for (i in 1:length(flabel)) 
      lines(cos(circle) * lab[i], sin(circle) * lab[i], lty=2, col="darkgrey")
    text(-lab/sqrt(2), lab/sqrt(2), cex=1.2, 
         labels = sprintf(paste0("%.", pmax(-digit, 0), "f"), label))
  }

  ## polygon set up
  xpoly = cos(new) * rep(df, each = m + 1)
  ypoly = sin(new) * rep(df, each = m + 1)
  circle2 = rev(circle)
  polygon(c(xpoly, cos(circle2) * radius), c(ypoly, sin(circle2) * radius), 
          col=col, border=border)
  segments(x0=cb*radius, y0=sb*radius, x1=x1, y1=y1, col=border)

  ## main title
  if(is.null(main)) {
    main = paste0(if(area.prop) "Area" else "Height", "-proportional ",
                  if(radius != 0) "Circular Histogram" else "Rose Diagram")
  }
  title(main=main)
}


  


