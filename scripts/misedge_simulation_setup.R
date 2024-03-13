library(ggplot2)
library(navigm)

N <- 200; P <- 100; Q <- 50; Q0 <- 2;
zeta <- qnorm(0.01); beta0 <- 0.4;

# 
dirname <-
  paste0(
    '~/simulation_rand_n_',
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
    '_misspecified_edge'
  )

#
dir.create(dirname)
setwd(dirname)
nrep <- 100

#
sig2_beta0 <- 0.1
sparsity_v <- c()

for(seed in 1:nrep){
  if(seed == 1){
    if (Q == Q0){
      nonzero_id <- 1:Q
    }else{
      nonzero_id <- sort(sample(1:Q, Q0))
    }
  }
  set.seed(seed)
  
  # set V 
  if(Q0 == 2){
    Vtmp <- rbind(matrix(rnorm(round(P * 0.3), -1, 1e-1), nrow = round(P * 0.3)),
               matrix(rnorm(round(P * 0.3), 0, 1e-1), nrow = round(P * 0.3)),
               matrix(rnorm((P-2 * round(P * 0.3)), 1, 1e-1), nrow = (P-2 * round(P * 0.3))))
    Vtmp <- cbind(Vtmp, rbind(matrix(rnorm(round(P * 0.3), 0, 1e-1), nrow = round(P * 0.3)),
                                matrix(rnorm(round(P * 0.3), 1, 1e-1), nrow = round(P * 0.3)),
                                matrix(rnorm((P-2 * round(P * 0.3)), 0, 1e-1), nrow = (P-2 * round(P * 0.3)))))
  }
  
  V <- matrix(rnorm(P*Q), nrow = P, ncol = Q)
  if(Q0==2){
    V[,nonzero_id] <- Vtmp
  }

  # plot V
  #
  tmp <- reshape::melt(V)
  ggplot(tmp, aes(X2, X1, fill=value))+
    geom_tile() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x='auxiliary variables', y='nodes')
  ggsave(paste0('V_heatmap',seed,'.pdf'), width = 4, height = 4)
  
  tmp <- dcast(tmp, X1~X2, value.var = "value")
  ggplot(tmp, aes(`1`,`2`)) +
    geom_point() +
    theme_bw() +
    labs(x = "variable 1", y = "variable 2")
  ggsave(paste0('V',seed,'.pdf'), width = 4, height = 4)
  
  #
  beta_true <- rep(0, Q)
  if(length(beta0)==1){
    if(beta0 > 0){
      beta_true[nonzero_id] <- rlnorm(Q0, log(beta0), sig2_beta0)
    }else{
      beta_true[nonzero_id] <- - rlnorm(Q0, log(-beta0), sig2_beta0)
    }
  }
    
  tmp <- Reduce ("+", lapply(1:Q, function(q){
    exp(-abs(matrix(V[,q], nrow = nrow(V), ncol = nrow(V)) 
        - t(matrix(V[,q], nrow = nrow(V), ncol = nrow(V))))) * beta_true[q]
  }))
  pe <- pnorm(zeta + tmp)
  
  pdf(paste0('ppi',seed,'.pdf'), height = 6, width = 4)
  hist(pe[upper.tri(pe)], xlab = 'edge inclusion probability', main = '')
  dev.off()
  
    delta <- matrix(rbinom(length(pe),1,pe), nrow = nrow(pe))
    delta[lower.tri(delta)] = t(delta)[lower.tri(delta)]
    diag(delta) <- 0
    # image(delta)
    sparsity_v <- c(sparsity_v, sum(delta[upper.tri(delta)])/sum(upper.tri(delta)))
    net <- simulate_data_from_adjacency_matrix(N,A=delta)
    
    
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

