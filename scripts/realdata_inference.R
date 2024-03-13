library(navigm)
#
option_list <- list(
  optparse::make_option(
    c("-s", "--seed"),
    type = "integer",
    default = NA,
    help = "seed for initialisation",
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
    default = '~/realdata/',
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
    c("-l", "--label"),
    type = "character",
    default = "",
    help = "data label",
    dest = "label"
  ),
  optparse::make_option(
    c("-f", "--fdr"),
    type = "numeric",
    default = 0.2,
    help = "false discovery rate for PPI filtering",
    dest = "fdr"
  )
)

#
args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))
list2env(args, envir = .GlobalEnv)

setwd(dir)

cat("Load ", paste0('realdata_cedar_fdr',fdr,
                     ifelse(label == "","",paste0("_",label)),'.rda'))
load(paste0('realdata_cedar_fdr',fdr,
            ifelse(label == "","",paste0("_",label)),'.rda'))

# 
lambda <- 2
v0_v <- s0_v <- seq(0.01, 1, length.out = 16)
v1 <- s1 <- 100
Q <- ncol(V)
P <- nrow(V)

#
list_hyper <- list(
  lambda = lambda,
  v0_v = v0_v,
  v1 = v1,
  a_tau = 2,
  b_tau = 2,
  a_sigma = 2,
  b_sigma = 2,
  a_o = 1,
  b_o = Q,
  a_rho = 1,
  b_rho = P)

#
if(is.na(seed)){
  
  # mu_beta is random
  #
  set.seed(123)
  list_init <-
    list(
      # mu_beta = rnorm(Q),
      mu_beta = rep(0, Q),
      sig2_inv_beta = rep(1, Q),
      alpha_tau = 1,
      beta_tau = 1,
      alpha_sigma = 1,
      beta_sigma = 1,
      alpha_o = 1,
      beta_o = Q-1,
      alpha_rho = 1,
      beta_rho = 1,
      tau1 = 1,
      tau2 = 1,
      o = 1/Q,
      rho = 1/P
    )
  list_init$beta <- list_init$mu_beta
  
}else{
  
  # random init for mu+_beta, alpha_tau, alpha_sigma
  set.seed(seed)
  
  list_init <-
    list(
      mu_beta = rnorm(Q),
      sig2_inv_beta = rep(1, Q),
      # provide enough flexibility to initialised means
      #
      alpha_tau = runif(1, 0.001, 2),
      beta_tau = 1,
      alpha_sigma = runif(1, 0.001, 2),
      beta_sigma = 1,
      #
      alpha_o = 1,
      beta_o = runif(1, 1 , Q-1),
      alpha_rho = 1,
      beta_rho = runif(1, 1 , P-1)
    )
  list_init$tau1 <- list_init$alpha_tau/list_init$beta_tau
  tau2 <- list_init$alpha_sigma/list_init$beta_sigma
  o <- list_init$alpha_o/list_init$beta_o
  rho <- list_init$alpha_rho/list_init$beta_rho
  
}

#
if(model !='GM'){
  version <- NULL
}

# inference
out <-
  navigm(Y = Y,
         V = V,
         method = model,
         inference = inference,
         list_hyper = list_hyper,
         list_init = list_init, 
         tol = 1e-3,
         maxit = 1e5,
         ne0 = c(mean_edge, stddev_edge^2),
         debug = T,
         full_output = T,
         version = version)

#
save(out,
     Y,
     V,
     file = paste0('out_',model,
                   ifelse(is.null(version),'',paste0('_v',version)),'_',
                   inference,ifelse(is.na(seed),'',paste0('_seed_', seed)),
                   '_zeta_mean_',mean_edge,
                   '_sd_',stddev_edge, 
                   ifelse(v0_v[1]!=1e4, paste0("_gridstart_",v0_v[1]), ""),
                   ifelse(label == "","",paste0("_",label)), 
                   ifelse(is.na(seed), "_init0",""),'.rda'))


