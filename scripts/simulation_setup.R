library(ggplot2)
library(navigm)
# Set up args for the simulation
#
# base
args <- list(corr= 0, col = NULL, blocks = NULL, dpar = NULL, codata = F, hub = F)
# sensitivity


# number of replicates
#
nrep <- 32
verbose <- T

# default values
#
if(! "N" %in% names(args)){
  args$N <- 200
}

if(! "P" %in% names(args)){
  args$P <- 100
}

if(! "Q" %in% names(args)){
  args$Q <- 50
}

if(! "Q0" %in% names(args)){
  args$Q0 <- 0
}

if(! "zeta" %in% names(args)){
  args$zeta <- -1.5
}

if(! "noise" %in% names(args)){
  args$noise <- 0.1
}

if(! "beta0" %in% names(args)){
  args$beta0 <- 0.5
}

if(! "corr" %in% names(args)){
  args$corr <- 0
}

list2env(args, envir = .GlobalEnv)

# if codata is included
# hub structure is on
if (codata) {
  
  if(hub == F){
    warning("hub structure exists if codata are included.")
    hub <- T
  }
  
}

# 
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
    zeta,
    '_noise_',
    noise,
    '_beta0_',
    beta0,
    '_codata_',
    codata,
    '_hub_',
    hub,
    ifelse(!is.null(args$corr),
           paste0('_corr_',corr),
           ""
    ),
    ifelse(!is.null(args$blocks),
           paste0('_blocks_',blocks),
           ""
    ),
    ifelse(!is.null(args$col),
           paste0('_col_',col),
           ""
    ),
    ifelse(!is.null(args$dpar),
           paste0('_l_', dpar),
           ""
    )
  )

#
dir.create(dirname)
setwd(dirname)


# correlation between codata
# 
if(!is.null(dpar)){
  
  # decreasing correlations with distance
  #
  # e.g. dpar = 1.5
  Sigma <- exp(-outer(1:Q,1:Q,'-')^2/(2*dpar^2))/2
  corr_index <- NULL
  # solve(Sigma)
  # image(Sigma)
  
}else if(!is.null(col)){
  
  # only several variables are correlated
  # 
  
  # e.g. col <- 5; corr <- 0.5
  corr_index <- sample(1:Q,col)
  Sigma <- diag(Q)
  Sigma[corr_index, corr_index] <-  corr
  diag(Sigma) <- 1
  # solve(Sigma)
  # image(Sigma)
  
}else if(!is.null(blocks)){
  
  # block correlation
  #
  # e.g. blocks <- 5; corr <- 0.5
  Sigma <- diag(Q)
  stopifnot( Q%%blocks == 0 )
  for(i in 1:(Q%/%blocks)){
    Sigma[((i-1) * blocks + 1) : (i * blocks),((i-1) * blocks + 1) : (i * blocks)] <- corr
  }
  diag(Sigma) <- 1
  corr_index <- NULL
  # solve(Sigma)
  # image(Sigma)
  
}else{
  
  # constant correlation
  #
  # e.g. corr <- 0.5
  Sigma <- matrix(corr, nrow = Q, ncol = Q)
  diag(Sigma) <- 1
  corr_index <- NULL
  # solve(Sigma)
  # image(Sigma)
  
}

#
sparsity_v <- c()
for(seed in 1:nrep){
  
  #
  set.seed(seed)
  
  # even codata = F generate V for algorithm testing. 
  
  # generate and plot codata
  #
  V <- generate_V(P, Q, min_gene = round(0.05 * P)) 
  if(verbose == T){
    print(quantile(apply(V, 2, function(x)sum(x>0.5))))
  }
 
  # plot V
  #
  tmp <- reshape::melt(V)
  ggplot(tmp, aes(X2, X1, fill=value))+
    geom_tile() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x='auxiliary variables', y='nodes')
  ggsave(paste0('V',seed,'.pdf'), width = 4, height = 4)
  
  pdf(paste0('Vhist',seed,'.pdf'),width = 6, height = 3)
  hist(cor(V)[upper.tri(cor(V))],xlim = c(-1,1), xlab = 'empirical correlations', main = 'correlation between auxiliary variables')
  dev.off()
  
  
  # fix nonzero_id across replicates
  #
  if(seed == 1){
    
    if (Q == Q0){
      
      nonzero_id <- 1:Q
      
    }else if(is.null(corr_index)){
      
      nonzero_id <- sort(sample(1:Q, Q0))
      
    }else{
      
      # only a few variables are correlated
      #
      # select 10% of correlated variables
      #
      ncorr_index <- max(1, round(length(corr_index) * 0.1))
      #
      if(Q0 > ncorr_index){
        nonzero_id <- sort(
          c(sample(corr_index, ncorr_index),
            sample(seq_len(Q)[-corr_index], Q0 - ncorr_index))
        )
        if(verbose){
          cat(ncorr_index, ' active variables are correlated and ', Q0 - ncorr_index, ' are uncorrelated.\n')
        }
      }else{
        nonzero_id <- sort(sample(corr_index, Q0))
      }
    }
  }
  
  

  if (codata) {
    
    # adjacency matrix
    #
    sig2_beta0 <- 0.1
    beta_true <- rep(0, Q)
    beta_true[nonzero_id] <- rlnorm(Q0, log(beta0), sig2_beta0)
    theta <- V %*% matrix(beta_true, ncol = 1)
    pe <- matrix(theta, nrow = P, ncol = P)
    pe <- pe + t(pe) + zeta
    pe <- pnorm(pe)
    delta <- 0 + (pe >= 0.5)
    diag(delta) <- 0
    
    # noise
    # 
    nnoise <- round(sum(delta[upper.tri(delta)]) * noise) 
    zero_entries <- which(upper.tri(delta) & delta==0, arr.ind = T)
    noise_id <- zero_entries[sample(nrow(zero_entries), nnoise),]
    noise_id <- matrix(noise_id, ncol=2)
    delta[noise_id] <- 1
    delta[noise_id[,c(2,1)]] <- 1
    
    # plot inclusion probs
    #
    pdf(paste0('ppi',seed,'.pdf'), height = 6, width = 4)
    hist(pe[upper.tri(pe)], xlab = 'edge inclusion probability', main = '')
    dev.off()
    
  } else{
    
    beta_true <- rep(0, Q)
    
    if (hub) {
      
      # random hub propensities
      #
      mu_theta <- 0.6
      v_theta <- 0.3
      theta <- rnorm(P, mu_theta, v_theta)
      pe <- matrix(theta, nrow = P, ncol = P)
      pe <- zeta + pe + t(pe)
      pe <- pnorm(pe)
      delta <- 0 + (pe >= 0.5)
      diag(delta) <- 0
      
    } else{
      
      # random networks
      #
      delta <- matrix(rbinom(P ^ 2, 1, pnorm(zeta)), nrow = P, ncol = P)
      diag(delta) <- 0
      delta[lower.tri(delta)] <- t(delta[upper.tri(delta)])
      
    }
  }
  
  # generate data
  # 
  net <- generate_data_from_adjancency(N,A=delta)
  
  
  # plot the simulated data
  #
  cat('Plot the simulated data ... \n')
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

  save(net, V, beta_true, file = paste0('data_seed', seed, '.rda'))
  
  bool_up <- upper.tri(net$A)
  sparsity_v <- c(sparsity_v, sum(net$A[bool_up] == 1)/sum(bool_up))
  
  if(verbose){
    cat('sparsity = ',sum(net$A[bool_up] == 1)/sum(bool_up), '\n')
  }
}


if(verbose){
  print(quantile(sparsity_v))
}

