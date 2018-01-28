MINDID_SFS <- function(P, scaleQ, j=2) {

  N   <- nrow(P)
  E_P <- ncol(P)

  sc_nbr <- length(scaleQ)
  Q_ni   <- vector("list",sc_nbr)
  Q_nbr  <- vector("numeric",sc_nbr)

  index    <- sc_nbr
  grp_cols <- names(P)
  for (nQ in rev(scaleQ)){
    r <- 1/nQ
    Q_ni[[index]] <- floor(P/r)[,list(count=.N),by=grp_cols]$count
    if (max(Q_ni[[index]])<= (j-1)) {
      stop('m is too large or there are not enough points')
    }
    Q_nbr[index] <- E_P*log(nQ)
    index <- index-1
  }

  logmMi <- vector('numeric',sc_nbr)
  for (i in 1:sc_nbr) {
    ni <- Q_ni[[i]]
    Q  <- Q_nbr[i]
    nMi <- 1
    NMi <- 1
    for(m in 1:(j-1)) {
      nMi <- nMi * (ni - m)
      NMi <- NMi * (N - m)
    }
    logmMi[i] <- Q*(j-1) + log(sum(ni * nMi) / (N * NMi))
  }

  delta <- log(sqrt((1/scaleQ)^2 * E_P))
  Mm    <- as.numeric(E_P-(-coef(lm(logmMi~delta))[2]/(j-1)))

  return(Mm)
}
