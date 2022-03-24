#' @title Multi-class Circular Density Curve
#' 
#' @description Function \code{cmdensity} can be used to plot 2-dimensional 
#'   density curves for circular data with multiple classes. The density curves 
#'   are stacked to avoid any overlap. 
#' 
#' @param funlist a list of functions which can be used to calculate the
#'   density values for each class, evaluated at given points defined by
#'   the first argument of the functions. The set of points is a sequence
#'   from \eqn{0} to \eqn{2\pi}, with length \code{n}.
#' @param funprop proportions for functions. It is 1 by default. A user can
#'   choose different proportions for the functions so as to represent
#'   different numbers of observations. If they do not add up to the number
#'   of functions (k), it will be normalised so that \code{sum(classprop) =
#'   k}.
#' @param radius the radius of the reference circle.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param total.area a positive number specifying the total area under all the
#'   density curves. If \code{total.area = NULL}, no scaling is applied, the
#'   plot is in the original scale. If \code{area.prop = TRUE}, the total area 
#'   is automatically unity without scaling.
#' @param n the number of points used to plot each density curve.  The
#'   larger the number is, the more accurate the curve is.
#' @param nlabels integer, for the number of levels to be plotted; if
#'   \code{0}, no label is plotted.
#' @param cols the colors to fill the area under each density curve, with
#'   the same order as the class.
#' @param borders the colors of the borders.
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
#' @param add logical; if \code{TRUE}, density curves are superimposed to
#'   the current plot, for example, the circular histograms, rose diagrams
#'   and stacked dot plots.
#' @param x.legend x coordinate to plot the legend.
#' @param y.legend y coordinate to plot the legend.
#' @param fill logical. If \code{TRUEt}, fills the regions with colors
#'   under/between the density curves. If \code{FALSE}, only the density
#'   curves are plotted.
#' @param lty line width
#' @param lwd line width
#' 
#' @concept stacked circular density curve
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
#' @seealso \code{\link{cdensity}}, \code{\link{cmhist}}
#' 
#' @importFrom graphics hist plot points text legend polygon lines title
#' @importFrom stats uniroot
#' @importFrom grDevices hcl
#' 
#' @export
#' 
#' @examples
#' # Load and pre-process the dataset
#' library(circular)
#' data("pigeons", package = "circular")
#' x = pigeons[,2] / 180 * pi                  # bearing
#' y = pigeons[,1]                             # treatment
#' vs = split(x, factor(y, unique(y)))    # list of classified value
#' prop = sapply(vs, length) / length(x)  # proportion of each class
#' 
#' # Define the kde function for each class using von Mises kernels
#' dvm = function(x, mu=0, kappa=1)  # von Mises density
#'   exp(kappa * cos(x - mu)) * (2 * pi * besselI(kappa, 0))^(-1)
#' kdevm = function(x, x0, bw=0.3) 
#'   rowMeans(outer(x, x0, dvm, 0.5 / (1 - exp(-bw^2 / 2))))
#' fs = list(function(x) kdevm(x, x0=vs[[1]]),
#'           function(x) kdevm(x, x0=vs[[2]]),
#'           function(x) kdevm(x, x0=vs[[3]]))
#' 
#' # stacked density curves for 3 classes
#' cmdensity(fs)                         # 1:1:1
#' cmdensity(fs, prop)                   # using proportions for functions
#'


cmdensity = function(funlist, funprop=1, radius=1/sqrt(base::pi),
                     area.prop=TRUE, total.area=1, 
                     n=500, nlabels=4, cols=NULL, borders=NULL, 
                     xlim=NULL, ylim=NULL, main=NULL,
                     type=c("null","compass","clock"), add=FALSE,
                     x.legend="bottomright", y.legend=NULL, fill=TRUE, lty=1,
                     lwd=1) {
  type = match.arg(type)
  pi = base::pi
  x = seq(0, 2 * pi, len = n)
  k = length(funlist)
  funprop = rep(funprop, len=k)
  funprop = funprop / sum(funprop)
  lty = rep(lty, len=k)
  lwd = rep(lwd, len=k)
  denmat = matrix(NA, nrow = n, ncol = k)
  for (i in 1:k) denmat[, i] = funlist[[i]](x)  # density values
  denmat = denmat * rep(funprop*k, rep(n, k))   # scaled by prop
  dencum = t(apply(denmat, 1, cumsum))
  denmax = dencum[, k]
  # factor = if(scale) scalefactor(denmax, radius, k) else 1
  if (is.null(total.area)) factor = 1 
  else factor = scalefactor(denmax, radius, total.area, area.prop) /
         (if(area.prop) k  else 1)
  dencumf = circtrans(dencum, radius, area.prop, factor)
  denmaxf = dencumf[,k]
  
  cx = cos(x)
  sx = sin(x)
  
  ## plot
  if (is.null(cols)) cols = hcl(seq(0, 240, len = k), c = 90, l = 80)
  cols = rep(cols, len=k)
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(cx * denmaxf))
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(sx * denmaxf))
  if(!add) {
    plot(0, type="n", asp=1, bty="n", axes=FALSE, ann=FALSE,
         xlim=xlim, ylim=ylim)
    points(0, 0, pch = 3)
    switch(type, 
           compass = text(cos(0:3*pi/2) * radius * 0.8, 
                          sin(0:3*pi/2) * radius * 0.8, 
                          labels = c("E", "N", "W", "S"),
                          cex = 1.2), 
           clock = {
             timeseq = seq(0, 2 * pi, len = 25)[-25]
             text(cos(timeseq) * radius * 0.85,  sin(timeseq) * radius * 0.85, 
                  labels = c(6:1, 24:7))
           },
           null =, )
  }
  
  ## labels
  if(!add && nlabels > 0) {
    ma = max(denmax) * 1.1
##    by = max(by.rounded(ma, nlabels), 0.01)
    by = signif(ma/nlabels, 1) 
    digit = floor(log10(by))
    flabel = seq(0, ma, by=by)[-1]             # area: density value for label
    lab = circtrans(flabel, radius, area.prop, factor)
    for (i in 1:length(flabel)) 
      lines(cos(x) * lab[i], sin(x) * lab[i], lty = 2, col="darkgrey")
    text(cex = 1.2, -lab/sqrt(2), lab/sqrt(2),
         labels = sprintf(paste0("%.", pmax(-digit, 0), "f"), flabel))
  }
  
  ## density curves
  for (i in 1:k) {
    den = dencumf[,k-i+1]
    if(fill) {
      polygon(cos(c(x, rev(x))) * c(den, rep(radius, n)),  # fill colors
              sin(c(x, rev(x))) * c(den, rep(radius, n)), 
              col=cols[k:1][i], border=cols[k:1][i])
      if(!is.null(borders))
        lines(cos(x) * den, sin(x) * den, lty=lty, lwd=lwd,
              col=borders[k:1][i])
    }
    else 
      lines(cos(x) * den, sin(x) * den, lty=lty, lwd=lwd, col=cols[k:1][i])
  }
  
  lines(cos(x) * radius, sin(x) * radius)  # reference circle
  if(!add) {
    ## main title
    if(is.null(main)) {
      main = paste0(if(area.prop) "Area" else "Height", "-proportional ",
                    "Stacked Circular Densities")
    }
    title(main = main)
    
    if(!is.null(x.legend)) {
      if(fill)
        legend(x.legend, y.legend, legend=paste0("f", 1:k), pch=22, pt.cex=2.2,
               pt.bg=cols, col={if(is.null(borders)) "black" else borders},
               xjust=0.5, yjust=0.5, bg="white")
      else 
        legend(x.legend, y.legend, legend=paste0("f", 1:k), lty=lty, lwd=lwd,
               col=cols, xjust=0.5, yjust=0.5, bg="white")
      }
  }
}






