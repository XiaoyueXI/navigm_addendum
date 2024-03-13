wlibrary(navigm)
#
option_list <- list(
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = 1,
    help = "simulation data seed (need to run simulation_setup first)",
    dest = "seed"
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
  ), 
  optparse::make_option(
    c("-a", "--v0"),
    type = "double",
    default = 0.142,
    help = "value of fixed v0",
    dest = "v0"
  ),
  optparse::make_option(
    c("-t", "--random_initial"),
    type = "logical",
    default = T,
    help = "random beta initialisation",
    dest = "bool_random_initial"
  )
)

#
args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))
list2env(args, envir = .GlobalEnv)

#
setwd(dir)
load(paste0('data_seed',seed,'.rda'))

# 
set.seed(seed)
lambda <- 2
s0 <- v0
v1 <- s1 <- 100
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
if(bool_random_initial){
  mu_beta <- rnorm(Q)
}else{
  mu_beta <- rep(0,Q)
}

list_init <-
  list(
    mu_beta = mu_beta,
    sig2_inv_beta = rep(1, Q),
    alpha_tau = 1,
    beta_tau = 1,
    alpha_sigma = 1,
    beta_sigma = 1,
    # -1 to make modes match EM
    alpha_o = 1,
    beta_o = Q-1,
    alpha_rho = 1,
    beta_rho = P-1,
    tau1 = 1,
    tau2 = 1,
    o = 1/Q,
    rho = 1/P
  )
list_init$beta <- list_init$mu_beta

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


save(out,
     net,
     V,
     beta_true,
     file = paste0('out_',model,ifelse(is.null(version),'',paste0('_v',version)),
                   '_',inference,
                   '_seed_', seed,
                   '_zeta_mean_',mean_edge, 
                   '_sd_',stddev_edge,
                   '_fixed_v0_',v0, 
                   ifelse(bool_random_initial, "","_init0"), '.rda'))


