#
library(navigm)
library(dplyr)
library(ggplot2)

#
dirname <- "~/uncertainty_simulation/"
dir.create(dirname)
setwd(dirname)

#
N <- 200; P <- 20; Q <- 10; Q0 <- 1
nonzero_id <- 5
zeta <- qnorm(0.01); beta0 <- 2

#
sparsity_v <- c()
net <- out_vbecm <- out_ecm <- V <- list()

#
for (seed in 1:100) {
  # simulate
  #
  V[[seed]] <- simulate_auxiliary_matrix(P, Q, empirical = F)
  beta_true <- rep(0, Q)
  beta_true[nonzero_id] <- beta0
  theta <- V[[seed]] %*% matrix(beta_true, ncol = 1)
  #
  pe <- matrix(theta, nrow = P, ncol = P)
  pe <- pe + t(pe) + zeta
  pe <- pnorm(pe)
  #
  delta <- matrix(rbinom(length(pe),1,pe), nrow = nrow(pe))
  delta[lower.tri(delta)] = t(delta)[lower.tri(delta)]
  diag(delta) <- 0
  #
  sparsity_v <- c(sparsity_v, print(sum(delta[upper.tri(delta)])/sum(upper.tri(delta))))
  net[[seed]] <- simulate_data_from_adjacency_matrix(N,A=delta)
  
  # inference
  #
  lambda <- 2
  grid_start <- 1e-2; grid_end <- 1
  v0_v <- s0_v <- seq(grid_start, grid_end, length.out = 16)
  v1 <- s1 <- 100
  
  #
  list_hyper <- list(
    lambda = lambda,
    v0_v = v0_v,
    v1 = v1,
    s0_v = s0_v, 
    s1 = s1,
    a_tau = 2,
    b_tau = 2,
    a_sigma = 2,
    b_sigma = 2,
    a_o = 1,
    b_o = Q,
    a_rho = 1,
    b_rho = P
   )
  
  #
  list_init <-
    list(
      # mu_beta = rnorm(Q),
      mu_beta = rep(0, Q),
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
  
  
  # inference
  #
  # vbecm
  #
  out_vbecm[[seed]] <-
    navigm(Y = net[[seed]]$Y,
           V = V[[seed]],
           method = "GMSS",
           inference = "VBECM",
           list_hyper = list_hyper,
           list_init = list_init, 
           tol = 1e-3,
           maxit = 1e5,
           # 
           ne0 = c(2, 2^2),
           debug = T,
           full_output = T,
           version = version)
  
  # ecm
  #
  out_ecm[[seed]] <-
    navigm(Y = net[[seed]]$Y,
           V = V[[seed]],
           method = "GMSS",
           inference = "ECM",
           list_hyper = list_hyper,
           list_init = list_init, 
           tol = 1e-3,
           maxit = 1e5,
           # 
           ne0 = c(2, 2^2),
           debug = T,
           full_output = T,
           version = version)
}

save(net, V, out_ecm, out_vbecm, file = "output.rda")


# load("output.rda")
#
# print(quantile(sparsity_v))
# zeta & beta
#
zeta_ecm <- lapply(out_ecm, function(x)x$estimates$zeta)
beta_ecm <- lapply(out_ecm, function(x)x$estimates$beta)
zeta_vbecm <- lapply(out_vbecm, function(x)x$estimates$mu_zeta)
sd_zeta_vbecm <- lapply(out_vbecm, function(x)sqrt(x$estimates$sig2_inv_zeta^(-1)))
beta_vbecm <- lapply(out_vbecm, function(x)x$estimates$mu_beta * x$estimates$m_gamma)
sd_beta_vbecm <- lapply(out_vbecm, function(x)sqrt(x$estimates$m_gamma * (x$estimates$sig2_inv_beta^(-1) + x$estimates$mu_beta^2) - x$estimates$m_gamma^2 * x$estimates$mu_beta^2))
Q <- length(beta_ecm[[1]])

# zeta
#
dfzeta <- rbind(data.frame(zeta = do.call("c",zeta_ecm), sd = 0, method = "ECM", id = 1:length(zeta_ecm)),
                data.frame(zeta = do.call("c",zeta_vbecm), sd = do.call("c",sd_zeta_vbecm), method = "VBECM", id = 1:length(zeta_vbecm)))

ggplot(dfzeta, aes(id, zeta, color=method)) + 
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = zeta - 1  * sd, ymax = zeta + 1  * sd), position = position_dodge(0.9))+
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "replicates", y = "estimated zeta", color = "") +
  geom_hline(yintercept = qnorm(0.01), linetype = 2)
ggsave("zeta_uq.pdf", width = 6, height = 4)

# 
dfzeta$inci <- F
dfzeta$inci[qnorm(0.01) >= dfzeta$zeta - 1.96  * dfzeta$sd & qnorm(0.01) <= dfzeta$zeta + 1.96  * dfzeta$sd] <- T
print(dfzeta %>% group_by(method) %>% summarise(sum(inci)))

# beta
#
dfbeta <- rbind(data.frame(beta = as.vector(do.call("rbind",beta_ecm)), 
                           sd = 0,
                           method = "ECM", 
                           id = rep(1:length(beta_ecm), times = Q),
                           rep = rep(1:Q, each = length(beta_ecm))),
                data.frame(beta = as.vector(do.call("rbind",beta_vbecm)), 
                           sd = as.vector(do.call("rbind",sd_beta_vbecm)),
                           method = "VBECM", 
                           id = rep(1:length(beta_vbecm), times = Q),
                           rep = rep(1:Q, each = length(beta_vbecm))))

nonzero_id <- 5; beta0 <- 2
beta_true <- rep(0,Q)
beta_true[nonzero_id] <- beta0

#
dfbetat <- data.frame(betat = beta_true, rep = 1:Q)
ggplot(dfbeta, aes(id, beta, color=method)) + 
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = beta - 1.96  * sd, ymax = beta + 1.96  * sd),position = position_dodge(0.9))+
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "replicates", y = "estimated beta", color = "") +
  geom_hline(data=dfbetat, aes(yintercept = betat), linetype = 2)+
  facet_wrap(rep~.,ncol = 2) 
ggsave("beta_uq.pdf", width = 12, height = 10)

#
dfbeta <- merge(dfbeta,dfbetat, by = "rep")
dfbeta$inci <- F
dfbeta$inci[dfbeta$betat >= dfbeta$beta - 1.96  * dfbeta$sd & dfbeta$betat <= dfbeta$beta + 1.96  * dfbeta$sd] <- T
print(dfbeta %>% group_by(method, rep) %>% summarise(sum(inci)))

