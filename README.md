# IDmining
Contains techniques for mining large high-dimensional data sets
by using the concept of Intrinsic Dimension (ID).

Version: 1.0.3

Author: Jean Golay and Mohamed Laib

Maintainer: Jean Golay Jean.Golay@unil.ch

URL: https://CRAN.R-project.org/package=IDmining

License: CC BY-NC-SA 4.0


# IDmining package installation: from github
install.packages("devtools")

devtools::install_github("jeangolay/IDmining")

library(IDmining)

# Example
bf <- Butterfly(1000)

require(colorRamps)
require(rgl)
c <- cut(bf$Y,breaks=64)
cols <- matlab.like(64)[as.numeric(c)]
plot3d(bf$X1,bf$X2,bf$Y,col=cols,radius=0.10,type="s",
xlab="",ylab="",zlab="",box=F)
axes3d(lwd=3,cex.axis=3)
grid3d(c("x+","y-","z"),col="black",lwd=1)
