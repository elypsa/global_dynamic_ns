est_global_cpp <- function(loc_f, glob_f, glob_param, loc_params) {

  states <- cbind(glob_f = glob_f, l1_glob_f = lag(glob_f),
                  l2_glob_f = lag(glob_f,2))
  rho_er <- loc_params['rho',]
  int_loc <- loc_params['beta0',]
  slp_loc <- loc_params['beta1',]

  K <- NCOL(loc_params)
  TT <- time(states[-1,])

  b <- matrix(int_loc * (1 - rho_er), nrow = K)
  H <- cbind(s1 = slp_loc, s2 = -slp_loc*rho_er)

  loc_f <- loc_level
  loc_f_star <- loc_f - sweep(lag(loc_f), 2, rho_er, FUN = "*")
  loc_f_star <- na.omit(loc_f_star)

  G <- matrix(0, 2,2)
  G[1,1] <- glob_param
  G[2,1] <- 1

  # Kalman filter
  # loc_f_star = b + H * states + e1 (R)
  # states = G*states + e2 (Q)

  P0 <- matrix(0, 2,2); P0[1,1] <- 1
  Q <- P0 # fixed for identification (see AR eq for glob fact)
  state_lag <- t(states[2, 1:2])
  P_lag <- P0
  R <- diag(loc_params['sigma2',]) # is sampled within chib alg


  ssm <- list(
    B0 = state_lag,
    P0 = P0,
    Dm = matrix(0, 2),
    Am = b,
    Fm = G,
    Hm = H,
    Qm = Q,
    Rm = R
  )

  kf <- kalman_filter(ssm, t(loc_f_star), smooth = T)

  # prepare containers:
  ss <- c(kf$B_tt[2,], kf$B_tt[1,NCOL(kf$B_tt)])
  states_carter <- xts(ss, order.by = time(glob_f))
  return(states_carter)
}



