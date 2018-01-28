#' Swiss Roll Data Set Generator
#'
#' Generates random points on the Swiss Roll manifold.
#' @usage SwissRoll(N=10000)
#' @param N The number of points to be generated (by default: \code{N = 10000}).
#' @return A \eqn{N \times 3}{N x 3} \code{data.frame} containing the
#' coordinates of the Swiss roll data points embedded in \eqn{\rm I\!R^3}{R^3}.
#' @references J. A. Lee and M. Verleysen (2007). Nonlinear Dimensionality Reduction, Springer, New York.
#' @examples
#' sim_dat <- SwissRoll(1000)
#' @importFrom stats runif
#' @export
SwissRoll <- function(N=10000) {

  if (N<=0) {stop('invalid argument')}

  x1 <- runif(N,min=-1,max=1)
  x2 <- runif(N,min=-1,max=1)
  x <- sqrt(2+2*x1)*cos(2*pi*sqrt(2+2*x1))
  y <- sqrt(2+2*x1)*sin(2*pi*sqrt(2+2*x1))
  z <- 2*x2

  dataset <- data.frame(x,y,z)

  return(dataset)
}
