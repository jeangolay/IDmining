#' Morisita-Based Filter for Regression Problems
#'
#' Executes the MBFR algorithm for supervised feature selection.
#'
#' @usage MBFR(XY, scaleQ, m=2, C=NULL)
#' @param XY A \eqn{N \times E}{N x E} \code{matrix}, \code{data.frame} or \code{data.table}
#' where \eqn{N} is the number of data points, \eqn{E} is the number of variables (i.e. the input
#' variables also called "features" + the output variable). The last column contains the
#' values of the output variable. And each variable (input + output) is rescaled to the
#' \eqn{[0,1]} interval by the function.
#' @param scaleQ A vector containing the values of \eqn{\ell^{-1}}{l^(-1)}
#' chosen by the user (see Details).
#' @param m The value of the parameter m (by default: \code{m=2}).
#' @param C The number of steps of the SFS procedure (by default: \code{C = E-1}).
#' @return A list of five elements:
#'  \enumerate{
#'  \item a vector containing the identifier numbers of the original features in the order
#'  they are selected through the Sequential Forward Selection (SFS) search procedure.
#'  \item the names of the corresponding features.
#'  \item the corresponding values of \eqn{Diss}.
#'  \item the ID estimate of the output variable.
#'  \item a \eqn{C \times 3}{C x 3} matrix containing: (column 1) the ID estimates of the subsets retained by the SFS
#'  procedure with the target variable; (column 2)  the ID estimates of the subsets retained by the
#'  SFS procedure without the output variable; (column 3) the values of \eqn{Diss} of the subsets
#'  retained by the SFS procedure.
#' }
#' @details
#'  \enumerate{
#'  \item \eqn{\ell}{l} is the edge length of the grid cells (or quadrats). Since the data
#'  (and consenquently the grid) are rescaled to the \eqn{[0,1]} interval, \eqn{\ell}{l} is equal
#'  to \eqn{1} for a grid consisting of only one cell.
#'  \item \eqn{\ell^{-1}}{l^(-1)} is the number of grid cells (or quadrats) along each axis of the
#'  Euclidean space in which the data points are embedded.
#'  \item \eqn{\ell^{-1}}{l^(-1)} is equal to \eqn{Q^{(1/E)}}{Q^(1/E)} where \eqn{Q} is the number
#'  of grid cells and \eqn{E} is the number of variables (or features).
#'  \item \eqn{\ell^{-1}}{l^(-1)} is directly related to \eqn{\delta}{delta} (see References).
#'  \item \eqn{\delta}{delta} is the diagonal length of the grid cells.
#'  \item The values of \eqn{\ell^{-1}}{l^(-1)} in \code{scaleQ} must be chosen according to the linear
#'  part of the \eqn{\log}{log}-\eqn{\log}{log} plot relating the \eqn{\log}{log} values of the
#'  multipoint Morisita index to the \eqn{\log}{log} values of \eqn{\delta}{delta} (or,
#'  equivalently, to the \eqn{\log}{log} values of \eqn{\ell^{-1}}{l^(-1)}) (see \code{logMINDEX}).
#' }
#' @author Jean Golay \email{Jean.Golay@@unil.ch}
#' @references
#' J. Golay, M. Leuenberger and M. Kanevski (2017). Feature selection for regression problems based
#' on the Morisita estimator of intrinsic dimension,
#' \href{http://www.sciencedirect.com/science/article/pii/S0031320317301905}{Pattern Recognition 70:126–138}.
#'
#' J. Golay, M. Leuenberger and M. Kanevski (2015).
#' \href{https://www.elen.ucl.ac.be/Proceedings/esann/esannpdf/es2015-41.pdf}{Morisita-based feature selection for
#' regression problems}.Proceedings of the 23rd European Symposium on Artificial Neural Networks, Computational
#' Intelligence and Machine Learning (ESANN), Bruges (Belgium).
#' @examples
#' \dontrun{
#' bf <- Butterfly(10000)
#'
#' fly_select <- MBFR(bf, 5:25)
#' var_order  <- fly_select[[2]]
#' var_perf   <- fly_select[[3]]
#'
#' dev.new(width=5, height=4)
#' plot(var_perf,type="b",pch=16,lwd=2,xaxt="n",xlab="",ylab="",
#'      ylim=c(0,1),col="red",panel.first={grid(lwd=1.5)})
#' axis(1,1:length(var_order),labels=var_order)
#' mtext(1,text = "Added Features (from left to right)",line = 2.5,cex=1)
#' mtext(2,text = "Estimated Dissimilarity",line = 2.5,cex=1)
#' }
#' @import data.table
#' @importFrom stats var lm coef
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MBFR <- function(XY, scaleQ, m=2, C=NULL) {

  if (!is.matrix(XY) & !is.data.frame(XY) & !is.data.table(XY)) {
    stop('XY must be a matrix, a data.frame or a data.table')
  }
  if (nrow(XY)<2){
    stop('at least two data points must be passed on to the function')
  }
  if (any(apply(XY, 2, var, na.rm=TRUE) == 0)) {
    stop('constant variables/features must be removed (they are not informative)')
  }
  if (!is.numeric(scaleQ) | length(scaleQ)<=1 | any(scaleQ<1) | any(scaleQ%%1!=0)) {
    stop('scaleQ must be a vector containing integers equal to or greater than 1')
  }
  if (length(m)!=1  | m<2 | m%%1!=0) {
    stop('m must be an integer equal to or greater than 2')
  }
  if (is.null(C)) {
     C <- ncol(XY)-1
  } else if (length(C)!=1  | C<1 | C > ncol(XY)-1 | C%%1!=0) {
     stop('C must be an integer between 1 and the number of input variables')
  }

  XY <- as.data.table(apply(XY, MARGIN = 2,
                            FUN = function(x) (x - min(x))/diff(range(x))))
  XY[XY==1] <- 1-0.5/max(scaleQ)

  E   <- ncol(XY)
  XsY <- XY[,E,with=F]
  X   <- XY[,1:(E-1),with=F]

  ID_Y <- round(MINDID_SFS(XsY, scaleQ, m),2)

  var_name_all <- colnames(X)
  var_order    <- character(C)
  mod_perf     <- matrix(NA,C,3)

  it <- txtProgressBar(min=0, max=C)
  for (i in 1:C){
    ID_temp <- matrix(NA,E-i,3)
    for (j in 1:(E-i)){
      xXsY <- cbind(X[,j,with=F],XsY)
      xXs  <- xXsY[,-ncol(xXsY),with=F]

      ID_temp[j,1] <- round(MINDID_SFS(xXsY, scaleQ, m),2)
      ID_temp[j,2] <- round(MINDID_SFS(xXs, scaleQ, m),2)

      ID_temp[j,3] <- round(ID_temp[j,1]-ID_temp[j,2],2)
    }
    diss_min     <- which.min(ID_temp[,3])
    var_order[i] <- colnames(X)[diss_min]
    mod_perf[i,] <- ID_temp[diss_min,]

    XsY <- cbind(X[,diss_min,with=F],XsY)
    X <- X[,-diss_min,with=F]

    setTxtProgressBar(it,i)
  }

  return(list(match(var_order,var_name_all),var_order,
              mod_perf[,3],ID_Y,mod_perf))
}
