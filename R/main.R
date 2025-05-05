##
rm(list=ls())
library(xts)
library(zoo)
library(reshape2)
library(dplyr)
library(kalmanfilter)
library(progress)

load('data/yield_curves.Rdata')
source('R/generate_latent_factors.R')
source('R/get_gl.R')
source('R/estimate_glob_ar1.R')
source('R/filter_global.R')

# estimate NS for each country
lambda <- 0.0609

loc_f <- lapply(yield_curves, function(x){
  res <- get_lf(yield_curve=x, lambda=lambda, curvature = F)
  return(res)
})


loc_level <- cbind(loc_f$US$level, loc_f$CA$level, loc_f$JP$level, loc_f$DE$level, loc_f$UK$level)
colnames(loc_level) <- c("US", "CA", "JP", "DE", "UK")
loc_slope <- cbind(loc_f$US$slope, loc_f$CA$slope, loc_f$JP$slope, loc_f$DE$slope, loc_f$UK$slope)
colnames(loc_slope) <- c("US", "CA", "JP", "DE", "UK")
# loc_curvature <- cbind(loc_f$US$curvature, loc_f$CA$curvature, loc_f$JP$curvature, loc_f$DE$curvature, loc_f$UK$curvature)
# colnames(loc_curvature) <- c("US", "CA", "JP", "DE", "UK")


# get global factors
# perform PCA
pr_level <- prcomp(loc_level, center = TRUE, scale. = TRUE)
pr_slope <- prcomp(loc_slope, center = TRUE, scale. = TRUE)
# pr_curvature <- prcomp(loc_curvature, center = TRUE, scale. = TRUE)

#extract global factors
glob_level <- as.xts(-pr_level$x[,"PC1"], order.by = time(loc_level))
glob_slope <- as.xts(-pr_slope$x[,"PC1"], order.by = time(loc_slope))
# glob_curvature <- as.xts(-pr_curvature$x[,"PC1"], order.by = time(loc_curvature))


reps <- 40e3
burn <- 20e3
pb <- progress_bar$new(total = reps, format = "  Progress [:bar] :percent in :elapsed")
for(i in 1:reps) {
  check_cond <- T
  counter <- 1
  while(check_cond) {
    counter <- counter + 1
    loc_params <- get_loc_params(loc = loc_level, global = glob_level,
                                 reps = 1, burn = 0)
    if(counter > 100) {
      loc_params['beta1','US'] <- -1*loc_params['beta1','US']
      print('counter > 100')
    }
    check_cond <- (loc_params['beta1', 'US'] <= 0)
  }

  global_param <- est_glob_ar1(glob_level)

  glob_level <- est_global_cpp(loc_level, glob_level, global_param, loc_params)

  pb$tick()

}


