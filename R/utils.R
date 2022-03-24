#################### some utility functions ####################

#' @title Circular Transformation Formula
#' 
#' @description The function performs circular transformation of density 
#'   or frequency, in an area-proportional or height-proportional manner.
#' 
#' @param x a numeric vector storing angular values between 0 and 2 pi, or
#'   an object that can be coerced to.
#' @param radius the radius of the reference circle.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' @param factor a positive number representing the scale factor to scale the
#'   entire plot.
#' 
#' @concept circular transformation
#' 
#' @return A numerical vector of the transformed values
#' 
#' @author Danli Xu <dxu452@aucklanduni.ac.nz>, Yong Wang <yongwang@auckland.ac.nz>
#' 
#' @references
#'
#' Xu, D. and Wang, Y. (2020) Area-proportional Visualization for
#' Circular Data. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{29}, 351-357.
#' 
#' @seealso \code{\link{scalefactor}}
#' 
#' @export
#' 
#' @examples
#' library(circular)
#' x = as.vector(rvonmises(20, circular(pi), 10))
#' circtrans(x)                            # area-proportional transformation
#' circtrans(x, area.prop = FALSE)         # height-proportional transformation
#' circtrans(x, factor = 2)                # with a scaling factor
#' 

circtrans = function(x, radius=0, area.prop=TRUE, factor=1) {
  if(area.prop) sqrt(2 * x * factor + radius^2)
  else x * factor + radius
}



#' @title Scaling Factor
#' 
#' @description The function calculates the scaling factor so that after scaling 
#'   the original density curve (before transformation), the total area after 
#'   transformation (excluding the reference circle) has the specified value.
#' 
#' @details Each value in x is a density value before transformation, for points
#'   equally-spaced on \eqn{[0,2\pi)}. For a smooth density curve, use a 
#'   reasonably large number of points, equally-spaced on \eqn{[0,2\pi)}. 
#'   The area under the density curve after transformation is then approximated 
#'   by that of the corresponding sectors. Note if \code{area.prop = TRUE}, 
#'   the scale factor is simply the value of \code{total.area}. 
#' 
#' @param x a numeric vector storing the heights of a density curve or
#'   a histogram.
#' @param radius the radius of the reference circle.
#' @param total.area a positive number specifying the total area.
#' @param area.prop logical; if \code{TRUE}, an area-proportional
#'   transformation is applied; if \code{FALSE}, a height-proportional
#'   transformationis applied.
#' 
#' @concept scale
#' 
#' @return A numerical value for the scaling factor
#' 
#' @author Danli Xu <dxu452@aucklanduni.ac.nz>, Yong Wang <yongwang@auckland.ac.nz>
#' 
#' @references
#'
#' Xu, D. and Wang, Y. (2020). Area-proportional Visualization for
#' Circular Data. \emph{Journal of Computational and Graphical
#' Statistics}, \bold{29}, 351-357.
#' 
#' @seealso \code{\link{circtrans}}
#' 
#' @export
#' 
#' @examples
#' dvm = function(x, mu=0, kappa=1)   # von Mises density
#' exp(kappa * cos(x - mu)) * (2 * pi * besselI(kappa, 0))^(-1)
#' x = dvm(seq(0, 2 * pi, len = 100), pi, 10)
#' 
#' scalefactor(x)                            # area-proportional transformation
#' scalefactor(x, area.prop = FALSE)         # height-proportional transformation
#' scalefactor(x, total.area = 2)            # total area of 2
#' scalefactor(x, area.prop = FALSE, total.area = 2)
#' 

## Finds the scaling factor so that after scaling the density curve (before
## transformation), the total area after transformation (excluding the
## reference circle) has the specified value.

## Each value in x is a density value before transformation, for points
## equally-spaced on [0,2*pi). For a smooth density curve, use a reasonably
## large number of points, equally-spaced on [0,2*pi). The area under the
## density curve after transformation is then approximated by that of the
## corresponding sectors.

## Each value in x is a density value before transformation, for points
## equally-spaced on \eqn{[0,2\pi)}. For a smooth density curve, use a reasonably
## large number of points, equally-spaced on [0,2*pi). The area under the
## density curve after transformation is then approximated by that of the
## corresponding sectors.

## x        vector, storing the heights of a density curve.


scalefactor = function(x, radius=0, total.area=1, area.prop=TRUE) {
  if(area.prop) return (total.area)
  xm = mean(x)
  x2m = mean(x^2)
  (sqrt(radius^2 * xm^2 + x2m / pi * total.area) - radius * xm) / x2m
}



## #########################
## 
## ## area of sectors (without the reference circle)
## ## x is a vector of density values evaluated at equally space points
## ## angle is the angle between two density values
## 
## sectorarea = function(x, rad=0, angle=pi) sum(x^2 - rad^2) / 2 * angle
## 
## 
## 
## ## Calculate (a part of) the area under the density curve.
## ## f is the density function.
## ## theta1 and theta2 lie inside the domain of f.
## ## If dimen = 1, the trapezoid rule is used.
## ## If dimen = 2, the area is the sum of sector areas. Each sector has a radius 
## ## equal to the density at the midpoint of that interval.
## ## n is the number of intervals to approximate the total area.
## 
## calarea = function(f, radius=0, area.prop=TRUE, factor=1, 
##                    theta1=0, theta2=2*pi, dimen=2, n=500) {
##   theta = seq(theta1, theta2, len = n)
##   mid = (theta[-n] + theta[-1]) / 2
##   val = circtrans(f(mid), radius, area.prop, factor)
##   if (dimen == 1) sum(f(theta[-n]) + f(theta[-1])) * (theta[2] - theta[1]) / 2
##   else sum(val^2 - radius^2) * (theta[2] - theta[1]) / 2
## }
## 
## ## example
## ## fcirc = function(x) rep(1, length(x))
## ## calarea(fcirc, 1)
## ## calarea(fcirc, 2)











## ## convert to area-proportional density values
## 
## areaprop = function(x, rad=0) sqrt(2 * x + rad^2)
## 
## 
## 
## ## halfstep = function(x) ceiling(x * 2) * 0.5    # round up to the nearst 0.5
## 
## ## height-prop only: scale the height of density curve so that total area equals to one
## fscale = function(k, x, rad=0) sum((x * k + rad)^2 - rad^2) * pi / length(x) - 1
## 
## ##
## 
## ## by.rounded = function(ma, nlevels) {
## ##   b = ma / nlevels
## ##   digits = ceiling(log10(b))
## ##   round(b, -(digits-1))
## ## }
## 
## 
## # cintegrate = function(f, lower=0, upper=2*pi, n=1001) {
## #   x = seq(lower, upper, len=n)
## #   dx = x[2] - x[1]       # angle
## #   x1 = x[-n] + dx / 2    # mid-points
## #   f1 = f(x1)
## #   sum(f1^2) * dx / 2
## # }
## #
## 
## ## height-prop only: scale the height of density curve 
## ## so that total area equals to a particular value
## fscaledata = function(k, x, rad=0, total=1) 
##   sum((x * k + rad)^2 - rad^2) * pi / length(x) - total
## 

