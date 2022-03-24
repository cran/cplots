#' @title Circular Stacked Dot Plot
#' 
#' @description  Function \code{cdotplot} can be used to plot 2-dimensional 
#'   stacked dot plot for circular data.
#' 
#' @param x a circular data object that is fully defined by the user.
#' @param nbins the number of bins of the circular histogram.  Internally,
#'   it is rounded to a multiple of 4.
#' @param radius the radius of the reference circle.  If \code{radius = 0},
#'   a rose diagram is produced; if \code{radius > 0}, a circular histogram
#'   is produced outside the reference circle.
#' @param unit the number of observations represented by each dot.  If
#'   \code{unit > 1}, it means that each dot represents multiple
#'   observations.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param total.area a positive number specifying the total area under the
#'   density curve. If \code{total.area = NULL}, no scaling is applied, the
#'   plot is in the original scale. If \code{area.prop = TRUE}, the total area 
#'   is automatically unity without scaling.
#' @param m the number of points within each bin to plot the circular dot
#'   plot. The larger the number is, the smoother the plot looks.
#' @param col the color to fill the bars.
#' @param border the color of the border around the bars.
#' @param xlim numeric vectors of length 2, giving the x coordinates
#'   ranges.
#' @param ylim numeric vectors of length 2, giving the y coordinates
#'   ranges.
#' @param main the main title (on top)
#' @param x.legend x coordinate to plot the legend.
#' @param y.legend y coordinate to plot the legend.
#' 
#' @details If the number of observations is relatively small, the
#'   usual circular stacked dot plot can be used with \code{unit = 1}.
#'   If the dataset is large, the dots may become too dense to
#'   visualize or count.  Setting \code{unit} to be any positive
#'   integer to allow each dot to represent more than one observation.
#'   If the number of observations in one bin is not a multiple of the
#'   specified unit, a partial dot can be used to represent the
#'   remainder at the top of the bin.
#' 
#' @concept circular stacked dot plot
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

#' @seealso \code{\link{cbarplot}}, \code{\link{cdensity}}, \code{\link{chist}}
#' 
#' @importFrom graphics hist plot points text polygon title
#' @importFrom stats uniroot
#' @importFrom circular rvonmises circular
#' 
#' @export
#' 
#' @examples
#' # 30 observations from two von Mises distributions
#' library(circular)
#' x = c(rvonmises(10, circular(pi/4), 5), rvonmises(20, circular(pi), 20))
#' cdotplot(x)                 # area-proportional dot plot
#' cdotplot(x, area = FALSE)   # height-proportional dot plot
#' 
#' # 900 observations from two von Mises distributions
#' y = c(rvonmises(300, circular(pi/4), 5), rvonmises(600, circular(pi), 20))
#' cdotplot(y, nbins=76, unit = 10)      # area-proportional (partial) dot plot 
#' cdotplot(y, nbins=76, unit = 10, area = FALSE) # height-proportional
#' 


cdotplot = function(x, nbins=36, radius=1, unit=NA, area.prop=TRUE, 
                    total.area=1, m=NA, col="lightblue", border="skyblue4", 
                    xlim=NULL, ylim=NULL, main=NULL,
                    x.legend="bottomright", y.legend=NULL) {
  x = as.vector(x)
  n = length(x)
  pi = base::pi
  nbins = max(4, round(nbins / 4) * 4)
  if(is.na(m)) m = max(ceiling(3600 / nbins), 2)
  br = seq(0, 2 * pi, len = nbins + 1)
  br2 = br[-(nbins+1)]
  cb2 = cos(br2)
  sb2 = sin(br2)
  hist = hist(x, breaks = br, plot = FALSE)
  if(is.na(unit)) unit = max(1, ceiling(n / 100))
  count = hist$counts / unit
  fr = hist$density    ## else fr = hist$counts
  # factor = if(scale) scalefactor(fr, radius) else 1
  factor = scalefactor(fr, radius, total.area, area.prop)
  if (is.null(total.area)) factor = 1 
  aa = diff(hist$breaks)[1] / 2  # length of semi-major
  bb = 0.5 / (2 * n * pi / nbins / unit)   # length of semi-minor
  he = m * 2     # no. points to draw half of ellipse circumference
  left = seq(pi / 2 * 3, pi / 2, len = he)
  dotarea = pi * aa * bb     # area of a full dot
  toparea = count %% 1 * dotarea        # area of the partial top dot
  
  # plot region
  fr2 = ceiling(fr * 2 * n * pi / nbins / unit) * unit / (2 * n * pi / nbins)
  fr2f = circtrans(fr2, radius, area.prop, factor)
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(cb2 * fr2f))
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(sb2 *  fr2f))
  plot(0, type="n", asp=1, axes=FALSE, ann=FALSE, xlim=xlim, ylim = ylim)
  
  # stacked dot
  for (i in 1:nbins) {
    ndot = ceiling(count[i])
    xc = hist$mids[i]
    yc = c(seq_len(ndot), rev(seq_len(ndot))) / 
      (2 * n * pi / nbins / unit) - bb
    si = rep(c(1, -1), rep(ndot, 2))
    xdot = rep(xc, he) + rep(si, rep(he, length(si))) * aa *
      cos(left)                                   # theta in 2d
    ydot = rep(yc, rep(he, length(yc))) + rep(si, rep(he, length(si))) * bb *
      sin(left)                                   # r in 2d
    tdot = c(xdot, seq(hist$mids[i], hist$mids[i] + 2 * pi / nbins, len = m))
    rdot = c(ydot, rep(0, m))
    frdot = circtrans(rdot, radius, area.prop, factor)
    polygon(cos(tdot) * frdot, sin(tdot) * frdot, col=col, border=border)
    if (count[i] %% 1 > 0) {
      toptheta = uniroot(farea, c(0 , 2 * pi), toparea[i], aa, bb)$root
      arc = seq(- (pi - toptheta) / 2, pi + (pi - toptheta) / 2, len = he * 2)
      tdot2 = rep(xc, he * 2) + cos(arc) * rtheta(arc, aa, bb)
      rdot2 = rep(yc[ndot], he * 2) + sin(arc) * rtheta(arc, aa, bb)
      tdot3 = seq(tdot2[length(tdot2)], tdot2[1], len = m)
      rdot3 = rep(rdot2[1], m)
      tdot23 = c(tdot2,tdot3)
      rdot23 = c(rdot2,rdot3)
      frdot23 = circtrans(rdot23, radius, area.prop, factor)
      polygon(cos(tdot23) * frdot23, sin(tdot23) * frdot23,
              col="white", border=border)
    }
  }

  angle = seq(0, 2*pi, len=500)      # draw reference circle
  polygon(cos(angle)*radius, sin(angle)*radius, col="white", border=border)
  points(0, 0, pch = 3)
  
  ## main title
  if(is.null(main)) 
    main = paste0(if(area.prop) "Area" else "Height", "-proportional ",
                  "Dot Plot")
  title(main=main)
  
  if(!is.null(x.legend)) {
    legend = substitute(unit == val, list(val=unit))
    legend(x.legend, y.legend, leg=legend, pch=21, cex=1.3, pt.bg=col,
           pt.cex=2.5, col=border, xjust=0.5, yjust=0.5)
  }
}


# length from a point on the circumference to the centre of an ellipse
rtheta = function(theta, major, minor) 
  major * minor / sqrt(sin(theta)^2 * major^2 + cos(theta)^2 * minor^2)

# area of a sector
fsector = function(theta, major, minor) 
  major * minor / 2 * (theta - atan((minor-major)*sin(2*theta) / 
                                      (minor+major+(minor-major)*cos(2*theta))))
# area of a segment of a sector
farea = function(t, s, major, minor) {
  t2 = (3 * pi + t) / 2
  t1 = (3 * pi - t ) / 2
  fsector(t2, major, minor) - fsector(t1, major, minor) - 
    rtheta(t2, major, minor) * rtheta(t1, major, minor) / 2 * sin(t) - s
}






