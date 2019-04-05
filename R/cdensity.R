#' @title Circular Density Curve
#' 
#' @description Function \code{cdensity} can be used to plot 2-dimensional 
#'   density curves for circular data.
#' 
#' @param f an R function that is to be plotted as a circular density or
#'   frequency.
#' @param radius the radius of the reference circle. If \code{radius = 0},
#'   no reference circle is produced, and the centre presents the point
#'   with zero density.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param total.area a positive number specifying the total area under the
#'   density curve. If \code{total.area = NULL}, no scaling is applied, the
#'   plot is in the original scale. If \code{area.prop = TRUE}, the total area 
#'   is automatically unity without scaling.
#' @param nlabels integer, for the number of levels to be plotted; if
#'   \code{0}, no label is plotted.
#' @param add logical; if \code{TRUE}, the density curve is superimposed to
#'   the current plot, for example, a circular histogram, a rose diagram or
#'   a stacked dot plot that has been produced in a similar manner.
#' @param n the number of points to plot the density curve.
#' @param col the color of the density line.
#' @param xlim numeric vectors of length 2, giving the x coordinates
#'   ranges.
#' @param ylim numeric vectors of length 2, giving the y coordinates
#'   ranges.
#' @param main the main title (on top)
#' 
#' @keywords circular density curve
#' 
#' @author Danli Xu <dxu452@aucklanduni.ac.nz>, Yong Wang <yongwang@auckland.ac.nz>
#' 
#' @references Xu, D. and Wang, Y. (2019) Area-proportional Visualization for 
#'   Circular Data (submitted).
#' 
#' @seealso \code{\link{cbarplot}}, \code{\link{cdotplot}}, \code{\link{chist}}
#' 
#' @importFrom graphics plot lines points title
#' @importFrom stats uniroot
#' @importFrom circular rvonmises circular
#' 
#' @export
#' 
#' @examples
#' # 600 observations from two von Mises distributions
#' library(circular)
#' x = c(rvonmises(200, circular(pi/4), 5), rvonmises(400, circular(pi), 20))
#' dvm = function(x, mu=0, kappa=1)   # von Mises density
#'   exp(kappa * cos(x - mu)) * (2 * pi * besselI(kappa, 0))^(-1)
#' f = function(x) 1/3 * dvm(x, pi/4, 5) + 2/3 * dvm(x, pi, 20)
#'
#' cdensity(f)                # plot the density in an area-proportional manner
#' 
#' chist(x)                   # circular histogram
#' cdensity(f, add=TRUE)      # superimpose the density curve
#' chist(x, area=FALSE)       # height-proportional circular histogram
#' cdensity(f, area=FALSE, add=TRUE)   # superimpose the density curve
#' 
#' chist(x, radius=0)                          # rose diagrams
#' cdensity(f, radius=0, add=TRUE)
#' chist(x, radius=0, area=FALSE)
#' cdensity(f, radius=0, area=FALSE, add=TRUE)
#' 


cdensity = function(f, radius=1/sqrt(base::pi), area.prop=TRUE, 
                    total.area=1, nlabels=4, add=FALSE, n=500, 
                    col="red", xlim=NULL, ylim=NULL, main=NULL) {
  pi = base::pi
  nlabels = max(0, ceiling(nlabels))
  x = seq(0, 2*pi, len=n)
  d = fval = f(x)
  # factor = if(scale) scalefactor(d, radius) else 1
  factor = scalefactor(d, radius, total.area, area.prop)
  if (is.null(total.area)) factor = 1 
  d = circtrans(d, radius, area.prop, factor)
  cx = cos(x)
  sx = sin(x)
  dx = cx * d
  dy = sx * d
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(dx))
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(dy))
  if(!add) {
    plot(0, type="n", xlim=xlim, ylim=ylim, asp=1, bty="n", axes=FALSE,
         ann=FALSE)
    if(nlabels > 0) {
      ma = max(fval) * 1.1
      by = signif(ma/nlabels, 1)
      flabel = seq(0, ma, by=by)[-1]           # area: density value for label
      lab = circtrans(flabel, radius, area.prop, factor)
      for (i in 1:length(flabel)) 
        lines(cx * lab[i], sx * lab[i], lty = 2, col="darkgrey")
      text(cex = 1.2, -lab/sqrt(2), lab/sqrt(2),
           labels = sprintf(if(by >= 0.1) "%.1f" else "%.2f", flabel))
    }
    lines(cx * radius, sx * radius, col="skyblue4")
    points(0, 0, pch = 3)
    if(is.null(main)) 
      main = paste0(if(area.prop) "Area" else "Height", "-proportional ",
                    "Circular Density Curve")
    title(main=main)
  }
  lines(dx, dy, col = col)
}





