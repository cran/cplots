#' @title Circular Bar Plot
#' 
#' @description Function \code{cbarplot} can be used to plot 2-dimensional 
#'   circular bar plots. The circular bar plots can only adopt the 
#'   height-proportional transformation because of the white space between bars.
#' 
#' @param x a numeric vector storing angular values between 0 and 2 pi, or
#'   an object that can be coerced to.
#' @param nbins the number of bins of the circular bar plot. Internally, it
#'   is rounded to a multiple of 4.
#' @param radius the radius of the reference circle.
#' @param prob logical; if \code{TRUE}, the circular histogram graphic is a
#'   representation of probability densities; if \code{FALSE}, a
#'   representation of frequencies.
#' @param nlabels integer, for the number of levels to be plotted; if
#'   \code{0}, no label is plotted
#' @param col the color to fill the bars.
#' @param border the color of the border around the bars.
#' @param m the number of points within each bin to plot the top of a
#'   bar. The larger the number is, the smoother the plot looks.
#' @param xlim numeric vectors of length 2, giving the x coordinates
#'   ranges.
#' @param ylim numeric vectors of length 2, giving the y coordinates
#'   ranges.
#' @param main the main title (on top)
#' 
#' @keywords circular bar plot
#' 
#' @author Danli Xu <dxu452@aucklanduni.ac.nz>, Yong Wang <yongwang@auckland.ac.nz>
#' 
#' @references Xu, D. and Wang, Y. (2019) Area-proportional Visualization for 
#'   Circular Data (submitted).
#' 
#' @seealso \code{\link{cdensity}}, \code{\link{cdotplot}}, \code{\link{chist}}
#' 
#' @importFrom graphics hist plot abline lines text points polygon title
#' @importFrom circular rvonmises circular
#' 
#' @export
#' 
#' @examples
#' # 600 observations from two von Mises distributions
#' library(circular)
#' x = c(rvonmises(200, circular(pi/4), 5), rvonmises(400, circular(pi), 20))
#'
#' cbarplot(x)                     
#' cbarplot(x, prob=FALSE)
#' cbarplot(x, radius=1, nlabels=0, col="lightblue")
#' cbarplot(x, radius=1, col="lightblue", border="skyblue4")
#' 


cbarplot = function(x, nbins=36, radius=1/sqrt(base::pi), prob=TRUE,
                    nlabels=4, col=NULL, border=NULL, m=NA, 
                    xlim=NULL, ylim=NULL, main=NULL) {
  x = as.vector(x)
  n = length(x)
  pi = base::pi
  nbins = max(4, round(nbins / 4) * 4)
  nlabels = max(0, ceiling(nlabels))
  if(is.na(m)) m = max(ceiling(360 / nbins), 2)
  br = seq(0, 2 * pi, len = nbins + 1)
  hist = hist(x, breaks = br, plot = FALSE)
  midangle = hist$mids    # same angle within each bar
  middens = hist$density  # same height within each bar
  halfbar = pi / nbins    # half intervel size of each bar = 2*pi / nbins / 2
  midseq = seq(-halfbar, halfbar, len = m)  # seq of the whole bar
  
  xmid = radius * cos(rep(midangle, each = m) + rep(midseq, nbins)) + 
    rep(middens, each = m) * cos(rep(midangle, each = m))
  ymid = radius * sin(rep(midangle, each = m) + rep(midseq, nbins)) + 
    rep(middens, each = m) * sin(rep(midangle, each = m))
  
  # add points on the circumference
  dim(xmid) = c(m, nbins)
  xco = rbind(xmid, radius * cos(br[-1]))
  dim(xco) = NULL
  dim(ymid) = c(m, nbins)
  yco = rbind(ymid, radius * sin(br[-1]))
  dim(yco) = NULL
  
  ## plot
  if(is.null(xlim)) xlim = c(-1,1) * max(abs(xco)) 
  if(is.null(ylim)) ylim = c(-1,1) * max(abs(yco)) 
  plot(0, type="n", asp=1, axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
  points(0, 0, pch = 3)
  if (nlabels > 0) {
    factor.fr = if(prob) 1 else 2 * pi * n / nbins
    ma = max(middens) * 1.1 * factor.fr
##    by = if(prob) max(by.rounded(ma, nlabels), 0.01)
#3         else max(by.rounded(ma, nlabels), 1)
    by = if(prob) signif(ma/nlabels, 1) else max(signif(ma/nlabels,1), 1)
    digit = floor(log10(by))
    label = seq(0, ma, by = by)[-1]
    hlabel = label / factor.fr + radius
    # abline(h=0, v=0, lty=2, col="darkgrey")
    circle = seq(0, 2*pi, len = 200)
    for (i in 1:length(hlabel))
      lines(cos(circle) * hlabel[i], sin(circle) * hlabel[i], lty=2,
            col="darkgrey")
    text(-hlabel/sqrt(2), hlabel/sqrt(2), cex = 1.2,
         labels=sprintf(paste0("%.", pmax(-digit, 0), "f"), label))
  }
  else {
    points(cos(pi * c(0.5, 1, 1.5, 2)) * radius, 
           sin(pi * c(0.5, 1, 1.5, 2)) * radius, pch = "+")
    text(cos(c(0, pi/2, pi, 3*pi/2)) * radius * c(0.85, 0.7), 
         sin(c(0, pi/2, pi, 3*pi/2)) * radius * c(0.85, 0.7), 
         cex = 1.3, labels = expression(0, frac(pi,2), pi, frac(3*pi,2)))
  }
  if(is.null(main)) main = "Circular Bar Plot"
  title(main=main)
  
  circle2 = seq(2 * pi, 0, len = 200)  
  polygon(c(xco, cos(circle2) * radius), c(yco, sin(circle2) * radius), 
          col=col, border=border)
  lines(cos(circle2) * radius,  sin(circle2) * radius)   # reference circle
}




