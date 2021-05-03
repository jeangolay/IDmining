#' The Multipoint Morisita Index for Spatial Patterns
#'
#' Computes the multipoint Morisita index for spatial patterns (i.e. 2-dimensional patterns).
#' @usage MINDEX_SP(X, scaleQ=1:5, mMin=2, mMax=5, Wlim_x=NULL, Wlim_y=NULL)
#' @param X A \eqn{N \times 2}{N x 2} \code{matrix}, \code{data.frame} or \code{data.table} containing the \eqn{X} and \eqn{Y}
#'   coordinates of \eqn{N} data points. The \eqn{X} coordinates must be given in the first column
#'   and the \eqn{Y} coordinates in the second column.
#' @param scaleQ  Either a single value or a vector. It contains the value(s) of \eqn{Q^{(1/2)}}{Q^(1/2)}
#'  chosen by the user where \eqn{Q} is the number of cells (or quadrats) of
#'  the \eqn{2D} grid (by default: \code{scaleQ = 1:5}).
#' @param mMin The minimum value of \eqn{m} (by default: \code{mMin = 2}).
#' @param mMax The maximum value of \eqn{m} (by default: \code{mMax = 5}).
#' @param Wlim_x A vector controlling the spatial extent of the \eqn{2D} gird along the
#'  \eqn{X} axis. It consists of two real values, i.e. \code{Wlim_x <- c(a,b)} where
#'  \eqn{b > a} (by default: \code{Wlim_x <- c(min(X[,1]),max(X[,1]))}).
#' @param Wlim_y A vector controlling the spatial extent of the \eqn{2D} gird along the
#'  \eqn{Y} axis. It consists of two real values, i.e. \code{Wlim_y <- c(a,b)} where
#'  \eqn{b > a} (by default: \code{Wlim_y <- c(min(X[,2]),max(X[,2]))}).
#' @return A \code{data.frame} containing the value of the m-Morisita index for each value of
#' \eqn{\delta}{delta} and \eqn{m}.
#' @details
#' \enumerate{
#' \item \eqn{Q^{(1/2)}}{Q^(1/2)} is the number of grid cells (or quadrats) along each of the two axes.
#' \item \eqn{Q^{(1/2)}}{Q^(1/2)} is directly related to \eqn{\delta}{delta} (see References).
#' \item \eqn{\delta}{delta} is the diagonal length of the grid cells.
#' }
#' @author Jean Golay \email{jeangolay@@gmail.com}
#' @references
#' J. Golay, M. Kanevski, C. D. Vega Orozco and M. Leuenberger (2014).
#' The multipoint Morisita index for the analysis of spatial patterns,
#' \href{https://www.sciencedirect.com/science/article/pii/S0378437114002659}{Physica A 406:191–202}.
#'
#' L. Telesca, J. Golay and M. Kanevski (2015). Morisita-based space-clustering analysis of Swiss seismicity,
#' \href{https://www.sciencedirect.com/science/article/pii/S0378437114008425}{Physica A 419:40–47}.
#'
#' L. Telesca, M. Lovallo, J. Golay and M. Kanevski (2016). Comparing seismicity declustering techniques
#' by means of the joint use of Allan Factor and Morisita index,
#' \href{https://link.springer.com/article/10.1007/s00477-015-1030-8}{Stochastic Environmental Research and Risk Assessment 30(1):77-90}.
#' @examples
#' sim_dat <- SwissRoll(1000)
#'
#' m <- 2
#' scaleQ <- 1:15 # It starts with a grid of 1^2 cell (or quadrat).
#'                # It ends with a grid of 15^2 cells (or quadrats).
#' mMI <- MINDEX_SP(sim_dat[,c(1,2)], scaleQ, m, 5)
#'
#' plot(mMI[,1],mMI[,2],pch=19,col="black",xlab="",ylab="")
#' title(xlab=expression(delta),cex.lab=1.5,line=2.5)
#' title(ylab=expression(I['2,'*delta]),cex.lab=1.5,line=2.5)
#'
#' \dontrun{
#' require(colorRamps)
#' colfunc <- colorRampPalette(c("blue","red"))
#' color <- colfunc(4)
#' dev.new(width=5,height=4)
#' plot(mMI[5:15,1],mMI[5:15,2],pch=19,col=color[1],xlab="",ylab="",
#'      ylim=c(1,max(mMI[,5])))
#' title(xlab=expression(delta),cex.lab=1.5,line=2.5)
#' title(ylab=expression(I['2,'*delta]),cex.lab=1.5,line=2.5)
#' for(i in 3:5){
#'   points(mMI[5:15,1],mMI[5:15,i],pch=19,col=color[i-1])
#' }
#' legend.text<-c("m=2","m=3","m=4","m=5")
#' legend.pch=c(19,19,19,19)
#' legend.lwd=c(NA,NA,NA,NA)
#' legend.col=c(color[1],color[2],color[3],color[4])
#' legend("topright",legend=legend.text,pch=legend.pch,lwd=legend.lwd,
#'        col=legend.col,ncol=1,text.col="black",cex=0.9,box.lwd=1,bg="white")
#'
#' xlim_l <- c(-5,5)     # By default, the spatial extent of the grid is set so
#' ylim_l <- c(-6,6)     # that it is the same as the spatial extent of the data.
#' xlim_s <- c(-0.6,0.2) # But it can be modified to cover either a larger (l)
#' ylim_s <- c(-1,0.5)   # or a smaller (s) study area (or validity domain).
#'
#' mMI_l <- MINDEX_SP(sim_dat[,c(1,2)], scaleQ, m, 5, xlim_l, ylim_l)
#' mMI_s <- MINDEX_SP(sim_dat[,c(1,2)], scaleQ, m, 5, xlim_s, ylim_s)
#' }
#' @import data.table
#' @importFrom stats var
#' @export
MINDEX_SP <- function(X, scaleQ=1:5, mMin=2, mMax=5, Wlim_x=NULL, Wlim_y=NULL){

  if (!is.matrix(X) & !is.data.frame(X) & !is.data.table(X) | ncol(X)!=2) {
    stop('X must be a matrix, a data.frame or a data.table with two columns')
  }
  if (ncol(X)!=2 | nrow(X)<2){
    stop('the X and Y coordinates must be passed on to the function (no more
          no less) with at least two data points')
  }
  if (any(apply(X, 2, var, na.rm=TRUE) == 0)) {
    stop('the X and/or Y coordinates must not be constant')
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
  if (is.null(Wlim_x)) {
     Wlim_x <- range(X[,1])
  } else if (!is.numeric(Wlim_x) | length(Wlim_x)!=2 | diff(Wlim_x)<=0){
     stop('the spatial extent of the grid is not set correctly.')
  }
  if (is.null(Wlim_y)) {
     Wlim_y <- range(X[,2])
  } else if (!is.numeric(Wlim_y) | length(Wlim_y)!=2 | diff(Wlim_y)<=0){
     stop('the spatial extent of the grid is not set correctly')
  }

  VD <- cbind(Wlim_x,Wlim_y) # set the extent of the grid

  id <- which(X[,1]>VD[2,1]| #
              X[,2]>VD[2,2]| # identification of the data that
              X[,1]<VD[1,1]| # do not fall into the grid
              X[,2]<VD[1,2]) #

  if (length(id)==nrow(X)){stop('no data points fall into the grid')}

  if (length(id)>0) {
    X_VD<-rbind(matrix(as.numeric(unlist(X[-id,])),
                        nrow=nrow(X[-id,])),VD)
  } else {
    X_VD<-rbind(matrix(as.numeric(unlist(X)),
                        nrow=nrow(X)),VD)
  }

  P_VD <- as.data.table(apply(X_VD, MARGIN = 2,
                        FUN = function(x) (x - min(x))/diff(range(x))))

  P  <- P_VD[-(nrow(P_VD):(nrow(P_VD)-1)),]

  P[P==1] <- 1-0.5/max(scaleQ) # replace 1 by a value sligthly below 1
  N <- nrow(P)                 # number of data points
  E <- ncol(P)                 # number of variables

  delta <- sqrt((diff(Wlim_x)/scaleQ)^2  # computation of the
               +(diff(Wlim_y)/scaleQ)^2) # values of delta

  sc_nbr <- length(scaleQ)
  Q_ni  <- vector("list",sc_nbr)
  Q_nbr <- vector("numeric",sc_nbr)

  index    <- sc_nbr
  grp_cols <- names(P)
  for (nQ in rev(scaleQ)){
    r <- 1/nQ
    Q_ni[[index]] <- floor(P/r)[,list(count=.N),by=grp_cols]$count
    if (max(Q_ni[[index]])<= (mMax-1)) {
      stop('mMax is too large or there are not enough points')
    }
    Q_nbr[index] <- nQ^E  # number of quadrats for each scale
    index <- index-1
  }

  mMi <- data.frame(Delta=delta)

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
      mMi[i,j-mMin+2] <- Q^(j-1) * sum(ni * nMi) / (N * NMi)
    }
    colnames(mMi)[j-mMin+2]<- paste("m",j,sep="")
  }

  return(mMi)
}
