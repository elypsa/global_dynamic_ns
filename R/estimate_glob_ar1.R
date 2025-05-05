est_glob_ar1 <- function(y) {
  x <- na.omit(dplyr::lag(y, 1))
  y_ <- y[2:length(y)]
  # freq version for simplicity

  B0 <- 0
  Sigma0 <- 1
  sigma2 <- 1 # set for identification

  B1 <- -5
  while(abs(B1)>=1) {
    V <- solve(solve(Sigma0) + (1/sigma2)*crossprod(x))
    M <- V%*%(solve(Sigma0)%*%B0+(1/sigma2)*crossprod(x, y_))
    B1 <- M + rnorm(1)*sqrt(V)
  }
  return(B1)


}
