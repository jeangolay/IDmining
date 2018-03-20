#' Functional Measure of Clustering Using the Morisita Estimator of ID
#'
#' Computes the functional m-Morisita index for a given set of threshold values.
#' @usage MINDID_FMC(XY, scaleQ, m=2, thd)
#' @param XY A \eqn{N \times E}{N x E} \code{matrix}, \code{data.frame} or \code{data.table}
#' where \eqn{N} is the number of data points and \eqn{E} is the number of variables
#' (i.e. the input variables + the variable measured at each measurement station). The last
#' column contains the variable measured at each measurement station. And each input variable
#' is rescaled to the [0,1] interval by the function. Typically, the input variables are the X
#' and Y coordinates of the measurement stations, but other or additional variables can be
#' considered as well.
#' @param scaleQ A vector containing the values of \eqn{\ell^{-1}}{l^(-1)}
#' chosen by the user (see Details).
#' @param m The value of the parameter m (by default: \code{m=2}).
#' @param thd Either a single value or a vector. It contains the value(s) of the threshold(s).
#' @return A \code{vector} containing the value(s) of the m-Morisita slope, \eqn{S_m}{Sm}, for each
#' threshold value.
#' @details
#' \enumerate{
#' \item \eqn{\ell}{l} is the edge length of the grid cells (or quadrats). Since the input variables
#' (and consenquently the grid) are rescaled to the \eqn{[0,1]} interval, \eqn{\ell}{l} is equal
#' to \eqn{1} for a grid consisting of only one cell.
#' \item \eqn{\ell^{-1}}{l^(-1)} is the number of grid cells (or quadrats) along each axis of the
#' Euclidean space in which the data points are embedded.
#' \item \eqn{\ell^{-1}}{l^(-1)} is equal to \eqn{Q^{(1/E)}}{Q^(1/E)} where \eqn{Q} is the number
#' of grid cells and \eqn{E} is the number of variables (or features).
#' \item \eqn{\ell^{-1}}{l^(-1)} is directly related to \eqn{\delta}{delta} (see References).
#' \item \eqn{\delta}{delta} is the diagonal length of the grid cells.
#' }
#' @author Jean Golay \email{Jean.Golay@@unil.ch}
#' @references
#' J. Golay, M. Kanevski, C. D. Vega Orozco and M. Leuenberger (2014).
#' The multipoint Morisita index for the analysis of spatial patterns,
#' \href{http://www.sciencedirect.com/science/article/pii/S0378437114002659}{Physica A 406:191–202}.
#'
#' J. Golay and M. Kanevski (2015). A new estimator of intrinsic dimension
#' based on the multipoint Morisita index,
#' \href{http://www.sciencedirect.com/science/article/pii/S0031320315002320}{Pattern Recognition 48 (12):4070–4081}.
#'
#' L. Telesca, J. Golay and M. Kanevski (2015). Morisita-based space-clustering analysis of Swiss seismicity,
#' \href{http://www.sciencedirect.com/science/article/pii/S0378437114008425}{Physica A 419:40–47}.
#' @examples
#' \dontrun{
#' bf    <- Butterfly(10000)
#' bf_SP <- bf[,c(1,2,9)]
#'
#' m      <- 2
#' scaleQ <- 5:25
#' thd    <- quantile(bf_SP$Y,probs=c(0,0.1,0.2,0.3,
#'                                    0.4,0.5,0.6,
#'                                    0.7,0.8,0.9))
#'
#' nbr_shuf    <- 100
#' Sm_thd_shuf <- matrix(0,length(thd),nbr_shuf)
#' for (i in 1:nbr_shuf){
#'   bf_SP_shuf      <- cbind(bf_SP[,1:2],sample(bf_SP$Y,length(bf_SP$Y)))
#'   Sm_thd_shuf[,i] <- MINDID_FMC(bf_SP_shuf, scaleQ, m, thd)
#' }
#' mean_shuf <- apply(Sm_thd_shuf,1,mean)
#'
#' dev.new(width=6, height=4)
#' matplot(1:10,Sm_thd_shuf,type="l",lty=1,col=rgb(1,0,0,0.25),
#'         ylim=c(-0.05,0.05),ylab=bquote(S[.(m)]),xaxt="n",
#'         xlab="",cex.lab=1.2)
#' axis(1,1:10,labels = FALSE)
#' text(1:10,par("usr")[3]-0.01,srt=45,ad=1,
#'      labels=c("0_100", "10_100","20_100","30_100",
#'               "40_100","50_100","60_100",
#'               "70_100","80_100","90_100"),xpd=T,font=2,cex=1)
#' mtext("Thresholds",side=1,line=3.5,cex=1.2)
#' lines(1:10,mean_shuf,type="b",col="blue",pch=19)
#'
#' legend.text<-c("Shuffled","mean")
#' legend.pch=c(NA,19)
#' legend.lwd=c(2,2)
#' legend.col=c("red","blue")
#' legend("topleft",legend=legend.text,pch=legend.pch,lwd=legend.lwd,
#'        col=legend.col,ncol=1,text.col="black",cex=1,box.lwd=1,bg="white")
#' }
#' @import data.table
#' @importFrom stats var lm coef
#' @export


MINDID_FMC <- function(XY, scaleQ, m=2, thd){

  if (!is.matrix(XY) & !is.data.frame(XY) & !is.data.table(XY)) {
    stop('XY must be a matrix, a data.frame or a data.table')
  }
  if (nrow(XY)<2){
    stop('at least two data points must be passed on to the function')
  }
  if (any(apply(XY, 2, var, na.rm=TRUE) == 0)) {
    stop('constant variables/features must be removed')
  }
  if (!is.numeric(scaleQ) | length(scaleQ)<=1 | any(scaleQ<1) | any(scaleQ%%1!=0)) {
    stop('scaleQ must be a vector containing integers equal to or greater than 1')
  }
  if (length(m)!=1 | m<2 | m%%1!=0) {
    stop('m must be an integer equal to or greater than 2')
  }
  if (!is.numeric(thd)) {
    stop('thd must be a scalar or a vector of real numbers')
  }

  XY_sc  <- as.data.table(cbind(apply(XY[,1:(ncol(XY)-1)],
                                      MARGIN = 2, FUN = function(XY)
                                        (XY - min(XY))/diff(range(XY))),
                                XY[,ncol(XY),drop=F]))

  Sm_thd <- vector('numeric',length(thd))
  for (i in 1:length(thd)){
    row_slct  <- as.logical(XY_sc[,ncol(XY_sc),with=F]>=thd[i])
    data_slct <- XY_sc[row_slct,1:(ncol(XY_sc)-1),with=F]
    Sm_thd[i] <- (ncol(XY_sc)-1)-MINDID_SFS(data_slct,scaleQ,m)
  }
  return(Sm_thd)
}
