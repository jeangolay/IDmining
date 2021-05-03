#' The Multipoint Morisita Index in 1, 2 or Higher Dimensions
#'
#' Computes the ln values of the multipoint Morisita index in 1, 2 or higher dimensional spaces.
#' @usage logMINDEX(X, scaleQ=1:5, mMin=2, mMax=2)
#' @param X A \eqn{N \times E}{N x E} \code{matrix}, \code{data.frame} or \code{data.table} where \eqn{N} is the number
#' of data points and \eqn{E} is the number of variables (or features). Each variable
#' is rescaled to the \eqn{[0,1]} interval by the function.
#' @param scaleQ  Either a single value or a vector. It contains the value(s) of \eqn{\ell^{-1}}{l^(-1)}
#' chosen by the user (by default: \code{scaleQ = 1:5}).
#' @param mMin The minimum value of \eqn{m} (by default: \code{mMin = 2}).
#' @param mMax The maximum value of \eqn{m} (by default: \code{mMax = 2}).
#' @return A \code{data.frame} containing the \eqn{\ln}{ln} value of the m-Morisita index for each value of
#' \eqn{\ln (\delta)}{ln(delta)} and \eqn{m}. Notice also that the values of
#' \eqn{\ln (\delta)}{ln(delta)} are provided with regard to the \eqn{[0,1]} interval.
#' @details
#' \enumerate{
#' \item \eqn{\ell}{l} is the edge length of the grid cells (or quadrats). Since the variables
#' (and consenquently the grid) are rescaled to the \eqn{[0,1]} interval, \eqn{\ell}{l} is equal
#' to \eqn{1} for a grid consisting of only one cell.
#' \item \eqn{\ell^{-1}}{l^(-1)} is the number of grid cells (or quadrats) along each axis of the
#' Euclidean space in which the data points are embedded.
#' \item \eqn{\ell^{-1}}{l^(-1)} is equal to \eqn{Q^{(1/E)}}{Q^(1/E)} where \eqn{Q} is the number
#' of grid cells and \eqn{E} is the number of variables (or features).
#' \item \eqn{\ell^{-1}}{l^(-1)} is directly related to \eqn{\delta}{delta} (see References).
#' \item \eqn{\delta}{delta} is the diagonal length of the grid cells.
#' }
#' @author Jean Golay \email{jeangolay@@gmail.com}
#' @references J. Golay and M. Kanevski (2015). A new estimator of intrinsic dimension
#' based on the multipoint Morisita index,
#' \href{https://www.sciencedirect.com/science/article/pii/S0031320315002320}{Pattern Recognition 48 (12):4070–4081}.
#' @examples
#' sim_dat <- SwissRoll(1000)
#'
#' m <- 2
#' scaleQ <- 1:15 # It starts with a grid of 1^E cell (or quadrat).
#'                # It ends with a grid of 15^E cells (or quadrats).
#' lnmMI <- logMINDEX(sim_dat, scaleQ, m, m)
#'
#' dev.new(width=5, height=4)
#' plot(exp(lnmMI[,1]),exp(lnmMI[,2]),pch=19,col="black",xlab="",ylab="")
#' title(xlab = expression(delta), cex.lab = 1.5,line = 2.5)
#' title(ylab = expression(I['2,'*delta]), cex.lab = 1.5,line = 2.5)
#'
#' dev.new(width=5, height=4)
#' plot(lnmMI[,1],lnmMI[,2],pch=19,col="black",xlab="",ylab="")
#' title(xlab = expression(paste("log(",delta,")")), cex.lab = 1.5,line = 2.5)
#' title(ylab = expression(paste("log(",I['2,'*delta],")")), cex.lab = 1.5,line = 2.5)
#' @import data.table
#' @importFrom stats var
#' @export
logMINDEX <- function(X, scaleQ=1:5, mMin=2, mMax=2) {

  if (!is.matrix(X) & !is.data.frame(X) & !is.data.table(X)) {
    stop('X must be a matrix, a data.frame or a data.table')
  }
  if (nrow(X)<2){
    stop('at least two data points must be passed on to the function')
  }
  if (any(apply(X, 2, var, na.rm=TRUE) == 0)) {
    stop('constant variables/features must be removed')
  }
  if (!is.numeric(scaleQ) | any(scaleQ<1) | any(scaleQ%%1!=0)) {
    stop('scaleQ must be an integer or a vector of integers equal to or
         greater than 1')
  }
  if (length(mMin)!=1 | length(mMax)!=1 | mMin<2 | mMax<2 | mMin%%1!=0 |
      mMax%%1!=0 | mMin>mMax) {
    stop('mMin and mMax must be integers equal to or greater than 2 and
          mMax must be equal to or greater than mMin')
  }

  P <- as.data.table(apply(X, MARGIN = 2,
                     FUN = function(x) (x - min(x))/diff(range(x))))
  N <- nrow(P)
  E <- ncol(P)

  delta <- sqrt((1/scaleQ)^2 * E)
  P[P==1] <- 1-0.5/max(scaleQ)

  sc_nbr <- length(scaleQ)
  Q_ni   <- vector("list",sc_nbr)
  Q_nbr  <- vector("numeric",sc_nbr)

  index    <- sc_nbr
  grp_cols <- names(P)
  for (nQ in rev(scaleQ)){
    r <- 1/nQ
    Q_ni[[index]] <- floor(P/r)[,list(count=.N),by=grp_cols]$count
    if (max(Q_ni[[index]])<= (mMax-1)) {
      stop('mMax is too large or there are not enough points')
    }
    Q_nbr[index] <- E*log(nQ)
    index <- index-1
  }

  logmMi <- data.frame(logDelta=log(delta))

  for (j in mMin:mMax) {
    for (i in 1:sc_nbr) {
      ni <- Q_ni[[i]]
      Q  <- Q_nbr[i]
      nMi <- 1
      NMi <- 1
      for(m in 1:(j-1)) {
        nMi <- nMi * (ni - m)
        NMi <- NMi * (N - m)
      }
      logmMi[i,j-mMin+2] <- Q*(j-1) + log(sum(ni * nMi) / (N * NMi))
    }
    colnames(logmMi)[j-mMin+2]<- paste("logm",j,sep="")
  }

  return(logmMi)
}
