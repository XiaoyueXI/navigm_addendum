library(navigm)
#
option_list <- list(
  optparse::make_option(
    c("-s", "--dseed"),
    type = "integer",
    default = 1,
    help = "simulation data seed (need to run simulation_setup first)",
    dest = "dseed"
  ),
  optparse::make_option(
    c("-t", "--iseed"),
    type = "integer",
    default = 1,
    help = "initialisation seed",
    dest = "iseed"
  ),
  optparse::make_option(
    c("-v", "--stddev"),
    type = "integer",
    default = 150,
    help = "priori knowledge on standard deviations of edge counts",
    dest = "stddev_edge"
  ),
  optparse::make_option(
    c("-m", "--mean"),
    type = "double",
    default = 50,
    help = "priori knowledge on average edge counts",
    dest = "mean_edge"
  ),
  optparse::make_option(
    c("-d", "--dir"),
    type = "character",
    default = '~',
    help = "directory that contains the simulated data",
    dest = "dir"
  ),
  optparse::make_option(
    c("-o", "--model"),
    type = "character",
    default = "GMSS",
    help = "the model to use (GM, GMN, GMSS)",
    dest = "model"
  ),
  optparse::make_option(
    c("-i", "--inference"),
    type = "character",
    default = "VBECM",
    help = "the inference method to use (VBECM, ECM)",
    dest = "inference"
  ),
  optparse::make_option(
    c("-n", "--version"),
    type = "integer",
    default = 2,
    help = "version of GM-ECM",
    dest = "version"
  )
)

#
args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))
list2env(args, envir = .GlobalEnv)

#
setwd(dir)
load(paste0('data_seed',dseed,'.rda'))

# 
lambda <- 2
v0 <- 0.142
v1 <- 100
s0 <- 1e-6
s1 <- 1
Q <- ncol(V)
P <- nrow(V)

#
list_hyper <- list(
  lambda = lambda,
  v0 = v0,
  v1 = v1,
  s0 = s0, 
  s1 = s1,
  a_tau = 2,
  b_tau = 2,
  a_sigma = 2,
  b_sigma = 2,
  a_o = 1,
  b_o = Q,
  a_rho = 1,
  b_rho = P)

#

#
set.seed(iseed)
list_init <-
  list(
    mu_beta = rnorm(Q),
    sig2_inv_beta = rgamma(1, 2, 2),
    alpha_tau = runif(1, 0.5, 5),
    beta_tau = runif(1, 0.5, 5),
    alpha_sigma = runif(1, 0.5, 5),
    beta_sigma = runif(1, 0.5, 5),
    # -1 to make modes match EM
    alpha_o = runif(1, 0.5, 5),
    beta_o = runif(1, 0.5, 2 * Q), # Q-1
    alpha_rho = runif(1, 0.5, 5),
    beta_rho = runif(1, 0.5, 2 * P)) # P-1


set.seed(iseed)
list_init <-
  list(
    mu_beta = rep(0,Q),
    sig2_inv_beta = rep(1,Q),
    alpha_tau = 1,
    beta_tau = 1,
    alpha_sigma = 1,
    beta_sigma = 1,
    # -1 to make modes match EM
    alpha_o = 1,
    beta_o = Q-1,
    alpha_rho = 1,
    beta_rho = P-1)


list_init$tau1 <- list_init$alpha_tau/list_init$beta_tau
list_init$tau2 <- list_init$alpha_sigma/list_init$beta_sigma
list_init$o = list_init$alpha_o/(list_init$alpha_o + list_init$beta_o)
list_init$rho = list_init$alpha_rho/(list_init$alpha_rho + list_init$beta_rho)
list_init$beta = list_init$mu_beta

#
if(model !='GM'){
  version <- NULL
}

# inference
out <-
  navigm_core(Y = net$Y,
              V = V,
              method = model,
              inference = inference,
              list_hyper = list_hyper,
              list_init = list_init, 
              tol = 1e-3,
              maxit = 1e5,
              ne0 = c(mean_edge, stddev_edge^2),
              debug = T,
              version = version)

#
save(out,
     net,
     V,
     beta_true,
     file = paste0('out_',model,ifelse(is.null(version),'',paste0('_v',version)),
                   '_',inference,
                   '_dseed_', dseed,
                   '_iseed_', iseed,
                   '_zeta_mean_',mean_edge, 
                   '_sd_',stddev_edge, 
                   '_v0_',v0,
                   '_v1_',v1,
                   '_s0_', s0, 
                   '_s1_', s1,'.rda'))


