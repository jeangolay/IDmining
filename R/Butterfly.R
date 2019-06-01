#' Butterfly Data Set Generator
#'
#' Generates a random simulation of the butterfly data set with a given number of points.
#' @usage Butterfly(N=10000)
#' @param N The number of points to be generated (by default: \code{N = 10000}).
#' @return A \eqn{N \times 9}{N x 9} \code{data.frame}. The first eight columns are the input variables,
#' and the last one is the output (or target) variable \eqn{Y}.
#' @author Jean Golay \email{jeangolay@@gmail.com}
#' @references J. Golay, M. Leuenberger and M. Kanevski (2016). Feature selection for regression problems based
#' on the Morisita estimator of intrinsic dimension,
#' \href{http://www.sciencedirect.com/science/article/pii/S0031320317301905}{Pattern Recognition 70:126–138}.
#' @examples
#' bf <- Butterfly(1000)
#'
#' \dontrun{
#' require(colorRamps)
#' require(rgl)
#'
#' c <- cut(bf$Y,breaks=64)
#' cols <- matlab.like(64)[as.numeric(c)]
#'
#' plot3d(bf$X1,bf$X2,bf$Y,col=cols,radius=0.10,type="s",
#'        xlab="",ylab="",zlab="",box=F)
#' axes3d(lwd=3,cex.axis=3)
#' grid3d(c("x+","y-","z"),col="black",lwd=1)
#'}
#' @importFrom stats runif
#' @export
Butterfly <- function(N=10000) {

  if (N<=0) {stop('invalid argument')}

  inpweight <- matrix(c(0.6655,1.2611,0.3961,-1.7065,0.8807,1.8260,1.3400,
                        1.2919,-1.3902,0.0743,0.8939,-0.3512,-1.7827,-0.5297,
                        1.9574,0.7962,1.5001,-0.4462,1.6856,1.5625),10,2)
  outweight <- c(1.3446,-0.0115,1.2770,0.5962,-0.8530,-0.7290,1.2339,0.1186,
                 0.5277,-0.6952)
  biashid   <- rep(0,10)

  cst_1<-0
  while(cst_1 == 0){
    x1      <- matrix(runif(N,-5,5),N,1)
    x2      <- matrix(runif(N,-5,5),N,1)
    x_core  <- cbind(x1,x2)
    if (min(x_core)>-5&max(x_core)<5){cst_1<-1}
  }

  nhid <- 10
  tmpHTest <- x_core %*% t(inpweight)
  biasMatrixTE <- matrix(rep(biashid,nrow(x_core)),nrow=nrow(x_core),
                         ncol=nhid,byrow=F)
  tmpHTest <- tmpHTest + biasMatrixTE
  HTest <- 1/(1 + exp(-1 * tmpHTest))
  Y <- HTest %*% outweight

  J3 <- log(x_core[,1]+5,10)
  J4 <- x_core[,1]^2-x_core[,2]^2
  J5 <- x_core[,1]^4+x_core[,2]^4
  # This is the correct version of J5 that is used in our research and papers
  # (i.e. an addition, not a subtraction).
  cst_2 <-0
  while(cst_2 == 0){
    I6 <- matrix(runif(N,-5,5),N,1)
    if (min(I6)>-5&max(I6)<5){cst_2<-1}
  }

  I7 <- log(I6+5,10)
  I8 <- I6+I7

  dataset <- data.frame(X1=x1,X2=x2,J3,J4,J5,I6,I7,I8,Y)

  return(dataset)
}
