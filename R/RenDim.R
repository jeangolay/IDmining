#' Renyi's Generalized Dimensions
#'
#' Estimates Rényi's generalized dimensions (or Rényi's dimensions of \eqn{qth}{qth} order). It is
#' mainly for \eqn{q=2}{q=2} that the result is used as an estimate of the intrinsic dimension of data.
#' @usage RenDim(X, scaleQ=1:5, qMin=2, qMax=2)
#' @param X A \eqn{N \times E}{N x E} \code{matrix}, \code{data.frame} or \code{data.table} where \eqn{N} is the number
#' of data points and \eqn{E} is the number of variables (or features). Each variable
#' is rescaled to the \eqn{[0,1]} interval by the function.
#' @param scaleQ A vector (at least two values). It contains the values of \eqn{\ell^{-1}}{l^(-1)}
#' chosen by the user (by default: \code{scaleQ = 1:5}).
#' @param qMin The minimum value of \eqn{q} (by default: \code{qMin = 2}).
#' @param qMax The maximum value of \eqn{q} (by default: \code{qMax = 2}).
#' @return A list of two elements:
#'  \enumerate{
#'   \item a \code{data.frame} containing the value of Rényi's information of \eqn{qth}{qth} order
#'   (computed using the natural logarithm) for each value of \eqn{\ln (\delta)}{ln(delta)}
#'   and \eqn{q}. The values of \eqn{\ln (\delta)}{ln(delta)} are provided with regard to the \eqn{[0,1]} interval.
#'   \item a \code{data.frame} containing the value of \eqn{D_q}{Dq} for each value of \eqn{q}.
#' }
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
#' @references
#' C. Traina Jr., A. J. M. Traina, L. Wu and C. Faloutsos (2000).
#' \href{https://periodicos.ufmg.br/index.php/jidm/article/view/4}{Fast feature
#' selection using fractal dimension}. Proceedings of the 15th Brazilian Symposium
#' on Databases (SBBD 2000), João Pessoa (Brazil).
#'
#' E. P. M. De Sousa, C. Traina Jr., A. J. M. Traina, L. Wu and C. Faloutsos (2007). A fast and
#' effective method to find correlations among attributes in databases,
#' \href{https://link.springer.com/article/10.1007/s10618-006-0056-4}{Data Mining and
#' Knowledge Discovery 14(3):367-407}.
#'
#' J. Golay and M. Kanevski (2015). A new estimator of intrinsic dimension
#' based on the multipoint Morisita index,
#' \href{https://www.sciencedirect.com/science/article/pii/S0031320315002320}{Pattern Recognition 48 (12):4070–4081}.
#'
#' H. Hentschel and I. Procaccia (1983). The infinite number of generalized
#' dimensions of fractals and strange attractors,
#' \href{https://www.sciencedirect.com/science/article/pii/016727898390235X}{Physica D 8(3):435-444}.
#' @examples
#' sim_dat <- SwissRoll(1000)
#'
#' scaleQ <- 1:15 # It starts with a grid of 1^E cell (or quadrat).
#'                # It ends with a grid of 15^E cells (or quadrats).
#' qRI_ID <- RenDim(sim_dat[,c(1,2)], scaleQ[5:15])
#'
#' print(paste("The ID estimate is equal to",round(qRI_ID[[1]][1,2],2)))
#' @import data.table
#' @importFrom stats var lm coef
#' @export
RenDim <- function(X, scaleQ=1:5, qMin=2, qMax=2) {

  if (!is.matrix(X) & !is.data.frame(X) & !is.data.table(X)) {
    stop('X must be a matrix, a data.frame or a data.table')
  }
  if (nrow(X)<2){
    stop('at least two data points must be passed on to the function')
  }
  if (any(apply(X, 2, var, na.rm=TRUE) == 0)) {
    stop('constant variables/features must be removed')
  }
  if (!is.numeric(scaleQ) | length(scaleQ)<=1) {
    stop('scaleQ must be a vector containing integers equal to or greater than 1')
  }
  if (length(qMin)!=1 | length(qMax)!=1) {
    stop('mMin and mMax must be real numbers')
  }

  P  <- as.data.table(apply(X, MARGIN = 2,
                      FUN = function(x) (x - min(x))/diff(range(x))))
  N <- nrow(P)
  E <- ncol(P)

  delta   <- sqrt((1/scaleQ)^2 * E)
  P[P==1] <- 1-0.5/max(scaleQ)

  sc_nbr <- length(scaleQ)
  Q_ni  <- vector("list",sc_nbr)

  index    <- sc_nbr
  grp_cols <- names(P)
  for (nQ in rev(scaleQ)){
    r <- 1/nQ
    Q_ni[[index]] <- floor(P/r)[,list(count=.N),by=grp_cols]$count
    index <- index-1
  }

  qRi <- data.frame(logDelta=log(delta))

  index <- 2
  for (j in qMin:qMax) {
    for (i in 1:sc_nbr) {
      ni  <- Q_ni[[i]]
      pri <- ni/N
      if (j == 1){
        qRi[i,index] <- -sum(pri * log(pri))
      } else {
        qRi[i,index] <- log(sum(pri^j)) * (1/(1 - j))
      }
    }
    colnames(qRi)[index] <- paste("RIq", j, sep="")
    index <- index+1
  }

  l  <- ncol(qRi)
  Dq <- vector("numeric",l-1)
  ID <- data.frame(qorder=qMin:qMax,Dq)

  for (i in 2:l){
    ID[i-1,2] <- -coef(lm(qRi[,i]~qRi[,1]))[2]
  }

  return(list(ID,qRi))
}
