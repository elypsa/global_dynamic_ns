# get global loadings

# get_gl_wrap <- function() {
#   list(level = get_gl(loc_level, glob_level),
#        slope = get_gl(loc_slope, glob_slope))
# }


get_loc_params <- function(loc, global, reps, burn) {
  apply(X = loc, MARGIN = 2, FUN = chib_alg, cbind(1, global), reps = reps, burn = burn)
}


chib_alg <- function(y, X, reps, burn) {
  n <- NCOL(X) # number of covariates
  B0 <- matrix(0, nrow = n, ncol = 1)
  Sigma0 <- diag(n) # sigma for B prior

  rho0 <- 0; rho <- rho0  # error autocorrelation coefficient
  Sigma0r <- 1 # sigma for rho

  T0 <- 0
  D0 <- 0
  sigma2 <- 1

  out <- matrix(NA, reps-burn, ncol = n + 2)
  colnames(out) <- c('rho', paste0('beta', 0:(n-1)), 'sigma2')


  for(i in 1:reps) {
    yStar <- na.omit(y - rho*lag(y,1))
    XStar <- na.omit(X - rho*lag(X,1))
    V <- solve(solve(Sigma0) + (1/sigma2)*crossprod(XStar))
    M <- V%*%(solve(Sigma0)%*%B0+(1/sigma2)*crossprod(XStar, yStar))

    B <- t(M) + rnorm(n)%*%chol(V)
    B <- t(B)

    y_eps <- y - X%*%B
    x_eps <- lag(y_eps, 1)

    y_eps <- y_eps[2:length(y_eps)]
    x_eps <- x_eps[2:length(x_eps)]

    rho_V <- solve(solve(Sigma0r) + (1/sigma2)*crossprod(x_eps))
    rho_M <- rho_V %*% (solve(Sigma0r)*rho0+(1/sigma2)*crossprod(x_eps,y_eps))


    rho <- 2
    while(rho >= 1) {
      rho <- c(rho_M + rnorm(1)*sqrt(rho_V))
    }

    resids <- yStar - XStar%*%B

    T1 <- T0 + nrow(XStar)
    D1 <- D0 + crossprod(resids)
    z0z0 <- crossprod(rnorm(T1))
    sigma2 <- c(D1/z0z0)

    if(i > burn) {
      out[i-burn,] <- c(rho, t(B), sigma2)
    }
  }
  return(colMeans(out))
}


test_chib <- function(size = 200, rho = 0.9, beta1 = 0.2, beta2 = 0.2, sigma2 = 0.5) {
  errors <- numeric(size)
  errors[1] <- rnorm(1, mean = 0, sd = sqrt(sigma2))

  for (i in 2:size) {
    errors[i] <- rho * errors[i - 1] + rnorm(1, mean = 0, sd = sqrt(sigma2))
  }
  X <- cbind(1, rnorm(size))
  y <- X %*% c(beta1, beta2) + errors

  out <- matrix(NA, 3, 4)
  rownames(out) <- c('DGP', 'OLS', 'Chib')
  colnames(out) <- c('rho', 'beta1', 'beta2', 'sigma2')
  out[1,] <- c(rho, beta1, beta2, sigma2)

  # ols
  ols_B <- solve(crossprod(X))%*%crossprod(X,y)
  out[2,] <- c(NA, ols_B[1], ols_B[2], NA)


  # chib
  chib_out <- chib_alg(y, X, reps = 20e3, burn = 10e3)
  out[3,] <- c(chib_out)
  return(out)
}
# test_chib(rho = 0)
# test_chib(size = 1e3, rho = 0.9)
