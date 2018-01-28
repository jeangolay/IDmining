#' Morisita-Based Filter for Redundancy Minimization
#'
#' Executes the MBRM algorithm for unsupervised feature selection.
#'
#' @usage MBRM(X, scaleQ, m=2, C=NULL, ID_tot=NULL)
#' @param X A \eqn{N \times E}{N x E} \code{matrix}, \code{data.frame} or \code{data.table}
#' where \eqn{N} is the number of data points and \eqn{E} is the number of variables (or features).
#' Each variable is rescaled to the \eqn{[0,1]} interval by the function.
#' @param scaleQ A vector containing the values of \eqn{\ell^{-1}}{l^(-1)}
#' chosen by the user (see Details).
#' @param m The value of the parameter m (by default: \code{m=2}).
#' @param C The number of steps of the SFS procedure (by default: \code{C = E}).
#' @param ID_tot The value of the full data ID if it is known a priori (by default:
#' the value of ID_tot is estimated using the Morisita estimator of ID witin
#' the function).
#' @return A list of four elements:
#'  \enumerate{
#'  \item a vector containing the identifier numbers of the original features in
#'  the order they are selected through the Sequential Forward Selection (SFS)
#'  search procedure.
#'  \item the names of the corresponding features.
#'  \item the corresponding ID estimates.
#'  \item the ID estimate of the full data set.
#'  }
#' @details
#'  \enumerate{
#'  \item \eqn{\ell}{l} is the edge length of the grid cells (or quadrats). Since the the variables
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
#' @examples
#' bf <- Butterfly(10000)
#'
#' bf_select <- MBRM(bf[,-9], 5:25)
#' var_order <- bf_select[[2]]
#' var_perf  <- bf_select[[3]]
#'
#' dev.new(width=5, height=4)
#' plot(var_perf,type="b",pch=16,lwd=2,xaxt="n",xlab="", ylab="",
#'      col="red",ylim=c(0,max(var_perf)),panel.first={grid(lwd=1.5)})
#' axis(1,1:length(var_order),labels=var_order)
#' mtext(1,text="Added Features (from left to right)",line=2.5,cex=1)
#' mtext(2,text="Estimated ID",line=2.5,cex=1)
#' @references
#' J. Golay and M. Kanevski (2017). Unsupervised feature selection based on the
#' Morisita estimator of intrinsic dimension,
#' \href{http://www.sciencedirect.com/science/article/pii/S0950705117303659}{Knowledge-Based Systems 135:125-134}.
#' @import data.table
#' @importFrom stats var lm coef
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
MBRM <- function(X, scaleQ, m=2, C=NULL, ID_tot=NULL) {

  if (!is.matrix(X) & !is.data.frame(X) & !is.data.table(X)) {
    stop('X must be a matrix, a data.frame or a data.table')
  }
  if (nrow(X)<2){
    stop('at least two data points must be passed on to the function')
  }
  if (any(apply(X, 2, var, na.rm=TRUE) == 0)) {
    stop('constant variables/features must be removed (they are not informative)')
  }
  if (!is.numeric(scaleQ) | length(scaleQ)<=1 | any(scaleQ<1) | any(scaleQ%%1!=0)) {
    stop('scaleQ must be a vector containing integers equal to or greater than 1')
  }
  if (length(m)!=1  | m<2 | m%%1!=0) {
    stop('m must be an integer equal to or greater than 2')
  }
  if (is.null(C)) {
    C <- ncol(X)
  } else if (length(C)!=1  | C<1 | C > ncol(X) | C%%1!=0) {
    stop('C must be an integer between 1 and ncol(X)')
  }

  P <- as.data.table(apply(X, MARGIN = 2,
                           FUN = function(x) (x - min(x))/diff(range(x))))
  P[P==1]      <- 1-0.5/max(scaleQ)
  E            <- ncol(P)
  var_name_all <- colnames(P)

  if (is.null(ID_tot)) {
    ID_tot <- round(MINDID_SFS(P, scaleQ, m),2)
  } else if (length(ID_tot)!=1 | ID_tot<0) {
    stop('ID_tot must be a non-negative real number')
  }

  var_name  <- colnames(P)
  var_order <- vector("character",C)
  mod_perf  <- vector("numeric",C)
  slct      <- c()
  tst       <- 1:ncol(P)

  it <- txtProgressBar(min=0, max=C)
  for (i in 1:C){
    ID_pi   <- vector("numeric",E+1-i)
    ID_diff <- vector("numeric",E+1-i)
    for (j in 1:(E+1-i)){
      temp       <- P[,c(slct,tst[j]),with=F]
      ID_pi[j]   <- round(MINDID_SFS(temp,scaleQ,m),2)
      ID_diff[j] <- round(abs(ID_tot-ID_pi[j]),2)
    }
    var_order[i] <- var_name[which.min(ID_diff)]
    mod_perf[i]  <- ID_pi[which.min(ID_diff)]

    slct     <- c(slct,tst[which.min(ID_diff)])
    tst      <- tst[-which.min(ID_diff)]
    var_name <- colnames(P)[tst]

    setTxtProgressBar(it,i)
  }

  return(list(match(var_order,var_name_all),var_order,mod_perf,ID_tot))
}
