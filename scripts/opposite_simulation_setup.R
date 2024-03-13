library(ggplot2)
library(navigm)
#
# cannot be too sparse; otherwise negative effects does not appear
#
args <- list(zeta = qnorm(0.2), beta0 = -0.5, Q0 = 3)
args <- list(zeta = qnorm(0.01),  beta0 = c(0.15, 0, -0.5), Q0 = 3)

#
list2env(args, envir = .GlobalEnv)

# set-up
#
N <- 200; P <- 100; Q <- 50; Q0 <- 3

dirname <-
  paste0(
    '~/simulation_n_',
    N,
    '_p_',
    P,
    '_q_',
    Q,
    '_q0_',
    Q0,
    '_zeta_',
    round(zeta,2),
    '_noise_',
    0,
    '_beta0_',
    ifelse(length(beta0) == 1, beta0, paste(beta0, collapse = "_")),
    '_opposite'
  )

#
dir.create(dirname)
setwd(dirname)

#
nrep <- 100
sparsity_v <- c()

for (seed in 1:nrep) {
  #
  set.seed(seed)
  
  #
  V <- simulate_auxiliary_matrix(P, Q, min_gene = round(0.05 * P)) 
  
  tmp <- reshape::melt(V)
  ggplot(tmp, aes(X2, X1, fill=value))+
    geom_tile() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x='auxiliary variables', y='nodes')
  ggsave(paste0('V_heatmap',seed,'.pdf'), width = 4, height = 4)
  
  #
  if(seed == 1){
    nonzero_id <- sort(sample(1:Q, Q0)) 
  }
  
  #
  sig2_beta0 <- 0.1
  beta_true <- rep(0, Q)
  if(length(beta0)==1){
    if(beta0 > 0){
      beta_true[nonzero_id] <- rlnorm(Q0, log(beta0), sig2_beta0)
    }else{
      beta_true[nonzero_id] <- - rlnorm(Q0, log(-beta0), sig2_beta0)
    }
  }else{
    if(length(beta0)!=Q0){
      stop("please provide beta0 as a number or a vector of length Q0")
    }else{
      beta_true[nonzero_id] <- beta0
    }
  }
  
  #
  theta <- V %*% matrix(beta_true, ncol = 1)
  pe <- matrix(theta, nrow = P, ncol = P)
  
  # 
  # if(all(beta0 < 0)){
  rzeta <- matrix(runif(length(pe),-(1-pnorm(zeta)), pnorm(zeta)), nrow = nrow(pe))
  # }else{
  #   rzeta <- matrix(runif(length(pe), - 2 * (1-pnorm(zeta)), 2 * pnorm(zeta)), nrow = nrow(pe))
  # }
  
  
  rzeta[lower.tri(rzeta)] = t(rzeta)[lower.tri(rzeta)]
  image(rzeta > 0)
  
  #
  pe <- pe + t(pe) + rzeta
  pe <- pnorm(pe)
  delta <- (pe >= 0.5)
  diag(delta) <- 0
  image(delta)
  
  #
  net <- simulate_data_from_adjacency_matrix(N,A=delta)
  sparsity_v <- c(sparsity_v, sum(net$A[upper.tri(net$A)])/sum(upper.tri(net$A)))
  
  # plot the simulated data
  #
  cat('Plot the simulated data ... \n')
  
  # adjacency matrix
  A <- reshape::melt(net$A)
  ggplot(A, aes(X1,X2, fill = value)) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = '', y = '', fill = 'A') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  ggsave(
    paste0('simulated_adjacency_matrix', seed, '.pdf'),
    width = 4.5,
    height = 4
  )
  
  # precision matrix
  Omega <- reshape::melt(net$Omega)
  ggplot(Omega, aes(X1, X2, fill = value)) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = '', y = '', fill = 'Omega') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  ggsave(
    paste0('simulated_precision_matrix', seed, '.pdf'),
    width = 4.5,
    height = 4
  )
  
  # empirical precision matrix if invertable
  try({
    EOmega <- reshape::melt(solve(cov(net$Y)))
    ggplot(EOmega, aes(X1,X2, fill = value)) +
      geom_tile() +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(x = '', y = '', fill = 'empirical \n Omega') +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())
    ggsave(
      paste0('simulated_empirical_precision_matrix', seed, '.pdf'),
      width = 4.5,
      height = 4
    )
  })
  
  # 
  save(net, V, beta_true, file = paste0('data_seed', seed, '.rda'))
  
}
quantile(sparsity_v)
print(mean(sparsity_v))

