library(navigm)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)

###############################

### low number of variables ###

###############################


dirname <- "~/simulation_n_200_p_100_q_3_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/"

setwd(dirname)

# all replicates
#
list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector('list',32)
list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector('list',32)

for(seed in 1:32){
  
  # gm
  #
  load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  bool_up <- upper.tri(net$A)
  list_GM[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GM[[seed]] <- net$A[bool_up]
  
  # gmn
  #
  load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMN[[seed]] <- net$A[bool_up]
  
  # gmss
  # 
  load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMSS[[seed]] <- net$A[bool_up]
  list_GMSS_c[[seed]] <- out$estimates$m_gamma
  list_true_GMSS_c[[seed]] <- beta_true!=0
  
}

# compare pROC curves using GM* and GMN
#

pdf('q3_proc_GM_GMN_GMSS.pdf',width = 4, height = 4)

par(pty = 's')

plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = 'black', main = 'pROC curves \n (edge selection)',lty = 3)
plot_roc(list_GMN, list_true_GMN, add = T, fpr_stop = 0.1, col = '#818181',lty = 2)
# plot_roc(list_GMSS, list_true_GMSS, add = T, fpr_stop = 0.1, col = '#D4D4D4',lty = 1)
legend("bottomright",     
       legend = c("GM*","GMN"),
       lty = c(3,2),
       cex = 0.7,
       col = c('black','#818181'))

dev.off()


# one replicate
#
seed <- 1
load(paste0('data_seed',seed,'.rda'))

df_true <- reshape::melt(net$A)
# df_margin_hub <- data.frame(margin1 = NA, margin2 = NA, X = 1:nrow(net$A))
# df_margin_hub$margin1 <- V %*% beta_true

#
load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
# set to NULL otherwise https://stackoverflow.com/questions/69666867/constant-warning-message-with-reshapemelt-in-r
rownames(out$estimates$m_delta) <- colnames(out$estimates$m_delta) <- NULL
df <- reshape::melt(out$estimates$m_delta)

#
load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
rownames(out$estimates$m_delta) <- colnames(out$estimates$m_delta) <- NULL
tmp <- reshape::melt(out$estimates$m_delta)

#
beta_est <- matrix(out$estimates$mu_beta, ncol=1)
df_margin_hub <- data.frame(margin = NA, X = 1:nrow(net$A))
df_margin_hub$margin <- V %*% beta_est

# gather GM & GMN edge estimation
#
df <- merge(df,tmp, by =c('X1','X2'))
df <- merge(df,df_true, by =c('X1','X2'))
df$GM <- df$GMN <- df$value.plot <- NA
df[df$value.x >= 0.5 & df$value == 1,]$GM <- 'TP'
df[df$value.x < 0.5 & df$value == 0,]$GM <- 'TN'
df[df$value.y >= 0.5 & df$value == 1,]$GMN <- 'TP'
df[df$value.y < 0.5 & df$value == 0,]$GMN <- 'TN'
df[is.na(df$GM),]$GM <- 'F'
df[is.na(df$GMN),]$GMN <- 'F'
df[df$X1 > df$X2,]$value.plot <- df[df$X1 > df$X2,]$GMN
df[df$X1 < df$X2,]$value.plot <- df[df$X1 < df$X2,]$GM
df[df$X1 == df$X2,]$value.plot <- NA

p <- ggplot(df, aes(x=X1, y=X2, fill=value.plot)) + 
  geom_tile() + 
  scale_fill_grey(start = 0.8, end = 0.2, na.value = 'white',na.translate = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(fill = '', x= '', y = '') +
  annotate("text", x= 25, y=75, label= "GM*", fontface="bold") + 
  annotate("text", x= 75, y=25, label= "GMN", fontface="bold") +
  theme(legend.position = 'top',
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm'))

p_legend <- ggpubr::get_legend(p)
p <- p + theme(legend.position = 'none')
px <- ggplot(df_margin_hub, aes(X, 1, fill = margin)) +
  geom_tile() +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_gradient(low = 'white', high = 'black')+
  labs(fill = 'hub propensity') +
  theme(legend.position="right",
        legend.direction = "horizontal",
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(color="white"),
        axis.ticks.x=element_line(color="white"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8))+
  guides(fill = guide_colourbar(barwidth = 5,  barheight = 0.3, title.position = 'top'))

px_legend <- ggpubr::get_legend(px)
px <- px + theme(legend.position = 'none')

#
df <- setDT(df)
df_margin_sum_tp <- df[X1 != X2, list(sum_GM = sum(GM == 'TP'), 
                                      sum_GMN = sum(GMN == 'TP')), by = 'X1']
df_margin_sum_tn <- df[X1 != X2, list(sum_GM = sum(GM == 'TN'), 
                                      sum_GMN = sum(GMN == 'TN')), by = 'X1']
df_margin_sum_tp  <- melt(df_margin_sum_tp, id.vars = 'X1')
df_margin_sum_tn  <- melt(df_margin_sum_tn, id.vars = 'X1')

py_tp <- 
  ggplot(df_margin_sum_tp, aes(X1, value, fill = factor(variable,c('sum_GM','sum_GMN'), c('GM*', 'GMN'))))+
  coord_flip() + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_grey() +
  theme_classic() + 
  labs(y = '# TP', fill ='') +
  theme(legend.position="top",
        legend.direction = "vertical",
        legend.box = "vertical",
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm'))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 

py_tp_legend <- ggpubr::get_legend(py_tp)
py_tp <- py_tp + theme(legend.position = 'none')

p + px + py_tp + p_legend + px_legend + py_tp_legend +
  plot_layout(design = "46
                        13
                        25",
              widths = c(0.8, 0.2),
              heights = c(0.05, 0.85, 0.05))

ggsave(paste0('q3_heatmap_GM_GMN_seed',seed,'.pdf'),width = 6, height = 5.5)

#################################################

### high number of variables  + null scenario ###

#################################################

# compare pROC curves using GM*, GMN and GMSS
#
for(label in c("q_50_q0_3","q_50_q0_0")){
  
  if (label == "q_50_q0_3"){
    
    dirname <- paste0("~/simulation_n_200_p_100_",label,"_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/")
    
  }else if(label == "q_50_q0_0"){
    
    dirname <- paste0("~/simulation_n_200_p_100_",label,"_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0/")
    
  }
  
  setwd(dirname)
  
  list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector('list',32)
  list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector('list',32)
  
  for(seed in 1:32){
    
    # gm
    #
    load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    bool_up <- upper.tri(net$A)
    list_GM[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GM[[seed]] <- net$A[bool_up]
    
    # gmn
    #
    load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMN[[seed]] <- net$A[bool_up]
    
    # gmss
    #
    load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
    
  }
  
  pdf(paste0(gsub('_','',label),'_proc_GM_GMN_GMSS.pdf'),width = 4, height = 4)
  
  par(pty = 's')
  
  plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = 'black', main = 'pROC curves \n (edge selection)',lty = 3)
  plot_roc(list_GMN, list_true_GMN, add = T, fpr_stop = 0.1, col = '#818181',lty = 2)
  plot_roc(list_GMSS, list_true_GMSS, add = T, fpr_stop = 0.1, col = '#D4D4D4',lty = 1)
  
  legend("bottomright",     
         legend = c("GM*","GMN","GMSS"),
         lty = 3:1,
         cex = 0.7,
         col = c('black','#818181','#D4D4D4'))
  
  dev.off()
  
  
  tmp <- do.call('rbind',list_GMSS_c)
  df_m_gamma <- data.table(M = apply(tmp, 2, mean),
                           SD = apply(tmp, 2, sd),
                           N = nrow(tmp),
                           ID = seq_len(ncol(tmp)))
  
  
  df_m_gamma[,L:=pmax(M-1.96*SD/sqrt(N), 0)]
  df_m_gamma[,U:=pmin(M+1.96*SD/sqrt(N),1)]
  tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = list_true_GMSS_c[[1]]!=0)
  df_m_gamma <- merge(df_m_gamma, tmp, by=c('ID'))
  
  ggplot(df_m_gamma, aes(x=ID, y=M, color=non_zero)) + 
    geom_point()+
    geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
    labs(x='\n  auxiliary variables', y = "posterior inclusion probability \n",
         color = 'active', title = 'Auxiliary variable PPIs\n')+
    theme_bw() +
    scale_x_discrete()+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))+
    scale_color_manual(breaks = c(FALSE,TRUE), 
                       labels = c('no','yes'),
                       values=c('#BABABA','black'))+
    theme(legend.position = 'bottom',panel.spacing = unit(1, "lines"))+
    geom_hline(aes(yintercept = 0.5), linetype = 2) +
    theme(plot.title =  element_text(hjust = 0.5, face = "bold"))
  
  ggsave(paste0(gsub('_','',label),'_PPI_GMSS_c.pdf'), width = 4, height = 3.5)
  
}

#####################

### null scenario ###

#####################

# estimated effects plot
#
# dirname <- '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
dirname <- '~/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0/'

setwd(dirname)
m1_beta_GMN <- m1_beta_GMSS <-  m1_zeta_GMN <- m1_zeta_GMSS <- list()
pauc_GMN <- pauc_GMSS <- c()

for(seed in 1:32){
  
  load(paste0('data_seed',seed,'.rda'))
  bool_up <- upper.tri(net$A)
  
  # gmn
  #
  load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  m1_beta_GMN[[seed]] <- out$estimates$mu_beta * out$estimates$m_gamma
  m1_zeta_GMN[[seed]] <- out$estimates$mu_zeta
  pauc_GMN <- c(pauc_GMN, compute_pauc(out$estimates$m_delta[bool_up], net$A[bool_up], fpr_stop = 0.1, standardise = T))
  
  # gmss
  #
  load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  m1_beta_GMSS[[seed]] <- out$estimates$mu_beta * out$estimates$m_gamma
  m1_zeta_GMSS[[seed]] <- out$estimates$mu_zeta
  pauc_GMSS <- c(pauc_GMSS, compute_pauc(out$estimates$m_delta[bool_up], net$A[bool_up], fpr_stop = 0.1, standardise = T))
  
}

# pauc
#
mean(pauc_GMN)
sd(pauc_GMN)/sqrt(32)
mean(pauc_GMSS)
sd(pauc_GMSS)/sqrt(32)

# plot effects
#
quantile(do.call('c',m1_zeta_GMN))
quantile(do.call('c',m1_zeta_GMSS))
m1_beta_GMN <- do.call('rbind', m1_beta_GMN)
m1_beta_GMSS <- do.call('rbind', m1_beta_GMSS)

df <- rbind(
  data.table(M = apply(m1_beta_GMN,2,mean),
             SD = apply(m1_beta_GMN,2,sd),
             N = apply(m1_beta_GMN,2,length),
             ID = seq_len(ncol(m1_beta_GMN)),
             method = 'GMN-VBECM'),
  data.table(M = apply(m1_beta_GMSS,2,mean),
             SD = apply(m1_beta_GMSS,2,sd),
             N = apply(m1_beta_GMSS,2,length),
             ID = seq_len(ncol(m1_beta_GMSS)),
             method = 'GMSS-VBECM')
)

df[,L:=M-1.96 * SD/sqrt(N)]
df[,U:=M+1.96 * SD/sqrt(N)]

dfs <- dcast(df, ID ~ method,value.var = c('M','L','U'))
dfs$active <- beta_true!=0
tmp <- range(c(df$M, df$L, df$U))

ggplot(dfs, aes(`M_GMSS-VBECM` , `M_GMN-VBECM`))+
  geom_errorbar(mapping = aes(ymin=`L_GMN-VBECM`, ymax=`U_GMN-VBECM`), color = 'grey') +
  geom_errorbarh(mapping = aes(xmin=`L_GMSS-VBECM`,xmax=`U_GMSS-VBECM`), color = 'grey') +
  geom_point() +
  theme_classic() +
  theme(legend.position = 'bottom') +
  theme(plot.title =  element_text(hjust = 0.5, face = "bold"))+
  labs(x='estimated effects GMSS', y='estimated effects GMN',
       title = 'Auxiliary variable effects\n') +
  coord_fixed()+
  scale_x_continuous(limits=tmp) +
  scale_y_continuous(limits=tmp) +
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0, linetype=2) 

ggsave('q50q00_scatter_effects_GMSS_GMN_c.pdf',width = 3, height = 3)

# # the previous run: seed 4 have substantially varying effects in GMN
# #
# load('~/simulation_n_200_p_100_q_50_q0_0_zeta_-1.84_noise_0_codata_FALSE_hub_FALSE_corr_0/out_GMN_VBEM_noanneal_select_AIC_seed_4_prior_zeta_prop0.01_stddev500.rda')
# tmp1 <- out$mu_beta
# tmp2 <- out$list_hyper
# tmp3 <- out$list_init
# # out$m_gamma
# out$mu_zeta
# 
# load('~/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0/out_GMN_VBECM_seed_4_zeta_mean_50_sd_150.rda')
# tmp1n <- out$estimates$mu_beta
# tmp2n <- out$args$list_hyper
# tmp3n <- out$args$list_init
# out$estimates$m_gamma
# out$estimates$mu_zeta


#########################

### inference methods ###

#########################

# GM
#
dirname <- '~/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0/'
setwd(dirname)

list_GM <- list_true_GM <- list_GM_ECM <- list_true_GM_ECM <- list()
pt_vbem <- index_vbem  <- it_vbem  <- vb_it_vbem <- list()
pt_em <- index_em  <- it_em  <- vb_it_em <- list()

for(seed in 1:32){
  
  # vbem
  #
  load(paste0('out_GM_v1_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  bool_up <- upper.tri(out$estimates$m_delta)
  list_GM[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GM[[seed]] <- net$A[bool_up]
  
  #
  tmp <- out$total_pt
  units(tmp) <- 'secs'
  pt_vbem[[seed]] <- tmp
  index_vbem[[seed]] <- out$index
  it_vbem[[seed]] <- sapply(out$full_output, function(x)x$it)
  vb_it_vbem[[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it)) # max iterations in vb
  
  # em
  #
  load(paste0('out_GM_v1_ECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  list_GM_ECM[[seed]] <- out$estimates$P1[bool_up]
  list_true_GM_ECM[[seed]] <- net$A[bool_up]
  
  #
  tmp <- out$total_pt
  units(tmp) <- 'secs'
  pt_em[[seed]] <- tmp
  index_em[[seed]] <- out$index
  it_em[[seed]] <- sapply(out$full_output, function(x)x$it)
  
}

pdf('proc_GM_GMSS_em_vbem.pdf',width = 8, height = 3)

par(mfrow = c(1,3), pty= 's')

plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = 'black',main='pROC curves \n (edge selection)',lty=2)
plot_roc(list_GM_ECM, list_true_GM_ECM, add = T, fpr_stop = 0.1, col = 'grey')

legend("bottomright",     
       legend = c("GM-VBECM","GM-ECM"),
       lty = c(2,1),
       cex = 0.7,
       col = c('black','grey'))

# runtime
#
round(mean(unlist(pt_em)),2)
round(sd(unlist(pt_em))/length(unlist(pt_em)),2)

round(mean(unlist(pt_vbem)),2)
round(sd(unlist(pt_vbem))/length(unlist(pt_vbem)),2)

# # why vbem is slower
# sapply(it_em, max)
# sapply(it_em, which.max)
# sapply(it_vbem, max)
# sapply(it_vbem, which.max)
# load(paste0('out_GM_v1_ECM_seed_',32,'_zeta_mean_50_sd_150.rda'))
# it_em[[32]]
# plot(out$full_output[[16]]$debugs$vec_ELBO_CM)
# load(paste0('out_GM_v1_VBECM_seed_',32,'_zeta_mean_50_sd_150.rda'))
# plot(out$full_output[[16]]$debugs$vec_ELBO_VBECM)
# plot(unlist(out$full_output[[16]]$debugs$list_ELBO))


# GMSS
#
dirname <- '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
setwd(dirname)

list_GMSS <- list_GMSS_c  <- list_true_GMSS <- list_true_GMSS_c <- vector('list',32)
list_GMSS_ECM <- list_GMSS_ECM_c <- list_true_GMSS_ECM <- list_true_GMSS_ECM_c <- vector('list',32)
pt_em <- index_em <- it_em <- vb_it_em <- list()
pt_vbem <- index_vbem <- it_vbem <- vb_it_vbem <- list()

for(seed in 1:32){
  
  # vbem
  #
  load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMSS[[seed]] <- net$A[bool_up]
  list_GMSS_c[[seed]] <- out$estimates$m_gamma
  list_true_GMSS_c[[seed]] <- beta_true!=0
  
  #
  tmp <- out$total_pt
  units(tmp) <- 'secs'
  pt_vbem[[seed]] <- tmp
  index_vbem[[seed]] <- out$index
  it_vbem[[seed]] <- sapply(out$full_output, function(x)x$it)
  vb_it_vbem[[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
  
  # em
  #
  load(paste0('out_GMSS_ECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  list_GMSS_ECM[[seed]] <- out$estimates$P1[bool_up]
  list_true_GMSS_ECM[[seed]] <- net$A[bool_up]
  list_GMSS_ECM_c[[seed]] <- out$estimates$P2
  list_true_GMSS_ECM_c[[seed]] <- beta_true!=0
  
  #
  tmp <- out$total_pt
  units(tmp) <- 'secs'
  pt_em[[seed]] <- tmp
  index_em[[seed]] <- out$index
  it_em[[seed]] <- sapply(out$full_output, function(x)sapply(x,function(y)y$it))
  
}

# edges
# 
plot_roc(list_GMSS, list_true_GMSS, add = F, fpr_stop = 0.1, col = 'black', 
         main = 'pROC curves \n (edge selection)', lty=2, ylab = '')
plot_roc(list_GMSS_ECM, list_true_GMSS_ECM, add = T, fpr_stop = 0.1, col = 'grey')
# abline(a=0,b=1)

legend("bottomright",
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2,1),
       cex = 0.7,
       col = c('black','grey'))


# auxiliary variables 
# 
plot_roc(list_GMSS_c, list_true_GMSS_c, add = F, fpr_stop = 0.1, col = 'black', 
         main = 'pROC curves \n (auxiliary variable selection)', ylab = '',lty = 2)
plot_roc(list_GMSS_ECM_c, list_true_GMSS_ECM_c, add = T, fpr_stop = 0.1, col = 'grey')

legend("bottomright",     
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2, 1),
       cex = 0.7,
       col = c('black','grey'))

dev.off()


# runtime
#
round(mean(unlist(pt_em)), 2)/60
round(sd(unlist(pt_em))/length(unlist(pt_em)), 2)
round(mean(unlist(pt_vbem)), 2)
round(sd(unlist(pt_vbem))/length(unlist(pt_vbem)), 2)

############################

### sensitivity analysis ###

############################

vec_dirnames <- c(
  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_200_p_50_q_50_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_100_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_100_p_50_q_50_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  
  '~/simulation_n_200_p_100_q_20_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_200_p_100_q_100_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  
  '~/simulation_n_200_p_100_q_50_q0_1_zeta_-0.75_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_200_p_100_q_50_q0_5_zeta_-2.1_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  
  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.2_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.3_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  
  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.75_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/',
  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.25_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
)

pauc_GMSS <- pauc_GMSS_c <-  pauc_GMN <-  pauc_GM <- list()

pt_GM <- index_GM <- it_GM <- vb_it_GM <- list()
pt_GMN <- index_GMN <- it_GMN <- vb_it_GMN  <- list() 
pt_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS  <- list()

sim_sparsity <- list()

for (dirname in vec_dirnames){
  
  setwd(dirname)
  
  list_GMSS <- list_true_GMSS <- list_GMN <- list_true_GMN <- list_GM <- list_true_GM <- list()
  list_GMSS_c <- list_true_GMSS_c <- list()
  
  pt_GM[[dirname]] <- index_GM[[dirname]]<- it_GM[[dirname]] <- vb_it_GM[[dirname]] <- list() 
  pt_GMN[[dirname]] <- index_GMN[[dirname]]<- it_GMN[[dirname]] <- vb_it_GMN[[dirname]] <- list()
  pt_GMSS[[dirname]] <- index_GMSS[[dirname]]<- it_GMSS[[dirname]] <- vb_it_GMSS[[dirname]] <- list()
  sim_sparsity[[dirname]] <-c()
  
  for(seed in 1:32){
    #
    load(paste0('data_seed',seed,'.rda'))
    bool_up <- upper.tri(net$A)
    sim_sparsity[[dirname]] <- c(sim_sparsity[[dirname]], sum(net$A[bool_up])/sum(bool_up))
    
    # gmss
    #
    load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    bool_up <- upper.tri(out$estimates$m_delta) 
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
    
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    pt_GMSS[[dirname]][[seed]] <- tmp
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMN[[seed]] <- net$A[bool_up]
    
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    pt_GMN[[dirname]][[seed]] <- tmp
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gm
    #
    load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    list_GM[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GM[[seed]] <- net$A[bool_up]
    
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    pt_GM[[dirname]][[seed]] <- tmp
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
  try({ pauc_GMSS[[dirname]] <- compute_pauc(list_GMSS,list_true_GMSS,  fpr_stop = 0.1, standardise = T) })
  try({ pauc_GMSS_c[[dirname]] <- compute_pauc(list_GMSS_c,list_true_GMSS_c,  fpr_stop = 0.1, standardise =T) })
  try({ pauc_GMN[[dirname]] <- compute_pauc(list_GMN,list_true_GMN,  fpr_stop = 0.1, standardise = T) })
  try({ pauc_GM[[dirname]] <- compute_pauc(list_GM,list_true_GM,  fpr_stop = 0.1, standardise = T) })
  
}

# sparsity
sapply(sim_sparsity, function(x)round(mean(x),4))

# check index and iterations
sapply(lapply(index_GM, unlist), quantile)
sapply(lapply(index_GM, unlist), function(x)sum(x==1))
sapply(lapply(index_GM, unlist), function(x)sum(x==16)) 
any(sapply(lapply(it_GM, unlist), max) == 1e5)
any(sapply(lapply(vb_it_GM, unlist), max)== 1e5)


sapply(lapply(index_GMN, unlist), quantile)
sapply(lapply(index_GMN, unlist), function(x)sum(x==1))
sapply(lapply(index_GMN, unlist), function(x)sum(x==16)) 
any(sapply(lapply(it_GMN, unlist), max) == 1e5)
any(sapply(lapply(vb_it_GMN, unlist), max)== 1e5)


sapply(lapply(index_GMSS, unlist), quantile)
sapply(lapply(index_GMSS, unlist), function(x)sum(x==1))
sapply(lapply(index_GMSS, unlist), function(x)sum(x==16)) 
any(sapply(lapply(it_GMSS, unlist), max) == 1e5)
any(sapply(lapply(vb_it_GMSS, unlist), max)== 1e5)

#
# create table for pauc
#
pauc_GMSS <- do.call('rbind',  pauc_GMSS)
pauc_GMSS <- apply(pauc_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GMSS_c <- do.call('rbind',  pauc_GMSS_c)
pauc_GMSS_c <- apply(pauc_GMSS_c, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GMN <- do.call('rbind',  pauc_GMN)
pauc_GMN <- apply(pauc_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GM <- do.call('rbind',  pauc_GM)
pauc_GM <- apply(pauc_GM, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(pauc_GM),t(pauc_GMN),t(pauc_GMSS), t(pauc_GMSS_c)), keep.rownames = T)
df$variable <- rep(c('pauc_GM','pauc_GMN','pauc_GMSS','pauc_GMSS_c'), each = length(unique(df$rn)))
setnames(df, c(paste0('V',1:3)), c('M','SD','N'))
df[,label:=paste0(round(M,2), ' (', round(SD/sqrt(N), 2)  ,') ')]
tmp <- dcast(df, rn~variable, value.var = 'label')
tmp[,rn := factor(rn, vec_dirnames)]

setkey(tmp, rn)
print(xtable(tmp[,c('rn','pauc_GM','pauc_GMN','pauc_GMSS','pauc_GMSS_c')]), include.rownames=FALSE)

#
# create table for runtime
#
pt_GMSS_v <- do.call('rbind',  lapply(pt_GMSS,'unlist'))
pt_GMSS_v <- apply(pt_GMSS_v, 1, function(x){c(mean(x), sd(x), length(x))})
pt_GMN_v <- do.call('rbind',  lapply(pt_GMN,'unlist'))
pt_GMN_v <- apply(pt_GMN_v, 1, function(x){c(mean(x), sd(x), length(x))})
pt_GM_v <- do.call('rbind',  lapply(pt_GM,'unlist'))
pt_GM_v <- apply(pt_GM_v, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(pt_GM_v),t(pt_GMN_v),t(pt_GMSS_v)), keep.rownames = T)
df$variable <- rep(c('time_GM','time_GMN','time_GMSS'), each = length(unique(df$rn)))
setnames(df, c(paste0('V',1:3)), c('M','SD','N'))
df[,label:=paste0(round(M,2), ' (', round(SD/sqrt(N), 2)  ,') ')]
tmp <- dcast(df, rn~variable, value.var = 'label')
tmp[,rn := factor(rn, vec_dirnames)]

setkey(tmp, rn)
print(xtable(tmp[,c('rn','time_GM','time_GMN','time_GMSS')]), include.rownames = F)


#######################

### runtime profile ###

#######################

# sample size

runtime_GM <- index_GM <- it_GM <- vb_it_GM <- list()
runtime_GMN <- index_GMN <- it_GMN <- vb_it_GMN <- list()
runtime_GMSS <-  index_GMSS <- it_GMSS <- vb_it_GMSS <- list()

for(n in c(20, 50, 100, 150, 200, 250, 300)){
  
  if(n!=100){
    dirname <-  paste0('~/simulation_n_',n,'_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0')
  }else{
    dirname <- '~/simulation_n_100_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0'
  }
  
  setwd(dirname)
  
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  for(seed in 1:32){
    
    # gm
    #
    load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]] <- c(index_GM[[dirname]] ,out$index)
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmss
    #
    load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

runtime_GM <- do.call('rbind', runtime_GM)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- do.call('rbind', runtime_GMN)
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMSS <- do.call('rbind', runtime_GMSS)
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})


df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c('GM','GMN','GMSS'), each = length(unique(df$rn)))
setnames(df, c(paste0('V',1:3)), c('M','SD','N'))
df[,L:=M-1.96 * SD/sqrt(N)]
df[,U:=M+1.96 * SD/sqrt(N)]
df[,N := as.numeric(str_extract(rn, "(?<=n_)\\d+"))]


ggplot(df, aes(x=N, y=M, group=variable, color = variable,lty=variable)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width =0.2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_y_continuous(trans='log10')+
  labs(x= 'sample size (N)', y='second (log-scale)',color ='',lty='')  +
  scale_color_grey()

ggsave('runtime_vs_n.pdf',width = 4, height = 3)


# graph size

runtime_GM <-   index_GM <- it_GM <- vb_it_GM <- list() 
runtime_GMN  <-index_GMN <- it_GMN <- vb_it_GMN <- list() 
runtime_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS <- list() 

for(p in seq(50, 300, by=50)){
  
  #
  if(p!=50 & p!=300){
    dirname <- paste0('~/simulation_n_200_p_',p,'_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0')
  }else{
    dirname <-paste0('~/simulation_n_200_p_',p,'_q_50_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0')
  }
  
  setwd(dirname)
  
  #
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  #
  for(seed in 1:32){
    # gm
    #
    load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn 
    #
    load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    
    # gmss
    #
    load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

runtime_GM <- do.call('rbind', runtime_GM)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- do.call('rbind', runtime_GMN)
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMSS <- do.call('rbind', runtime_GMSS)
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c('GM','GMN','GMSS'), each = length(unique(df$rn)))
setnames(df, c(paste0('V',1:3)), c('M','SD','N'))
df[,L:=M-1.96 * SD/sqrt(N)]
df[,U:=M+1.96 * SD/sqrt(N)]
df[,P := as.numeric(str_extract(rn, "(?<=p_)\\d+"))]

ggplot(df, aes(x=P, y=M, group=variable, color = variable,lty=variable)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width =0.2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_y_continuous(trans='log10')+
  labs(x= 'graph size (P)', y='second (log-scale)',color ='',lty='')  +
  scale_color_grey()
ggsave('runtime_vs_p.pdf',width = 4, height = 3)


# breakdown
#
dirname <- '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
setwd(dirname)
pt <-list()

for(seed in 1:32){
  #
  load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  pt[[seed]] <- 
    sapply(out$full_output, function(x){
      tmp <- x$pt
      units(tmp) <- 'secs'
      return(tmp)
    })
}

#
pt <- do.call('cbind',pt)
pt <- apply(pt, 1, function(x)c(mean(x),sd(x),length(x)))
pt <- data.table(t(pt))
pt$v0 <- seq(1e-4,1,length.out = 16)
colnames(pt)[1:3] <- c('M','SD','N')
pt[,L:= M-1.96 * SD/sqrt(N)]
pt[,U:= M+1.96 * SD/sqrt(N)]

ggplot(pt, aes(v0, M)) +
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U),width = 0.02) +
  theme_bw() +
  scale_y_continuous(trans = 'log10') + 
  labs(x='\n spike standard deviation', y = 'second (log-scale)') +
  theme(
    plot.margin = margin(5.5,5.5,25.5,5.5)
  )
ggsave('runtime_GMSS_reference.pdf',width = 4,height = 3)


# number of auxiliary variables
#
runtime_GM <-   index_GM <- it_GM <- vb_it_GM <- list() 
runtime_GMN  <-index_GMN <- it_GMN <- vb_it_GMN <- list() 
runtime_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS <- list() 

for(q in c(20, 50, 100)){
  
  #
  if(q!=50){
    dirname <- paste0('~/simulation_n_200_p_100_q_',q,'_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0')
  }else{
    dirname <-paste0('~/simulation_n_200_p_100_q_',q,'_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0')
  }
  
  setwd(dirname)
  
  #
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  for(seed in 1:32){
    # gm
    #
    load(paste0('out_GM_v2_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    load(paste0('out_GMN_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmss
    #
    load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
    tmp <- out$total_pt
    units(tmp) <- 'secs'
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

runtime_GM <- do.call('rbind', runtime_GM)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- do.call('rbind', runtime_GMN)
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMSS <- do.call('rbind', runtime_GMSS)
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})


df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c('GM','GMN','GMSS'), each = length(unique(df$rn)))
setnames(df, c(paste0('V',1:3)), c('M','SD','N'))
df[,L:=M-1.96 * SD/sqrt(N)]
df[,U:=M+1.96 * SD/sqrt(N)]
df[,Q := as.numeric(str_extract(rn, "(?<=q_)\\d+"))]

ggplot(df, aes(x=Q, y=M, group=variable, color = variable,lty=variable)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width =0.2) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_y_continuous(trans='log10')+
  scale_x_continuous(breaks=c(20,50,100)) + 
  labs(x= 'number of variables (Q)', y='second (log-scale)',color ='',lty='')  +
  scale_color_grey()

ggsave('runtime_vs_q.pdf',width = 4, height = 3)


########################

### model selection ###

########################

dirname <- '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
setwd(dirname)

df_model_criteria <- data.table()

for (seed in 1:32) {
  load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_50_sd_150.rda'))
  bool_up <- upper.tri(out$estimates$Omega)
  Y <- scale(net$Y, center = TRUE, scale = FALSE) # center the data
  S <- crossprod(Y)
  N <- nrow(Y)
  
  #
  estimates <-lapply(out$full_output, function(x){x$estimates})
  aic <- sapply(estimates, function(x) AIC_GSS(x, N))
  bic <- sapply(estimates, function(x) BIC_GSS(x, N))
  ebic <- sapply(estimates, function(x) EBIC_GSS(x, N = N))
  pauc <- sapply(out$full_output, function(x){
    compute_pauc(x$estimates$m_delta[bool_up], net$A[bool_up],fpr_stop = 0.1, standardise = T)
  })
  paucc <- sapply(out$full_output, function(x){
    compute_pauc(x$estimates$m_gamma, beta_true!=0,fpr_stop = 0.1, standardise = T)
  })
  
  #
  tmp <- data.table(aic = aic, 
                    bic = bic, 
                    ebic = ebic,
                    pauc = pauc, 
                    paucc = paucc,
                    v0 = out$arg$list_hyper$v0_v,
                    seed = seed)
  
  df_model_criteria <- rbind(df_model_criteria, tmp)
  
}

tmp <- df_model_criteria[,list(aic = which.min(aic),
                               pauc = which.max(pauc),
                               paucc = which.max(paucc)
                               ),by = 'seed']

df_model_criteria <- melt(df_model_criteria, id.vars = c('v0', 'seed'))
df_model_criteria_summary <- df_model_criteria[,list(M = mean(value),
                                                     SD = sd(value),
                                                     N = length(value)), by = c('v0', 'variable')]
df_model_criteria_summary[,L:= M - 1.96 * SD/sqrt(N)]
df_model_criteria_summary[,U:= M + 1.96 * SD/sqrt(N)]
df_model_criteria_summary[, minM := min(M), by = 'variable']
df_model_criteria_summary[, maxM := max(M), by = 'variable']
df_model_criteria_summary[, type := 'general']
df_model_criteria_summary[!grepl('pauc',variable) & M == minM, type := 'best']
df_model_criteria_summary[variable == 'pauc' & M == maxM, type := 'best']

#
ggplot(df_model_criteria_summary[!grepl('auc',variable),], aes(x = round(v0,3), y = M, color = type)) + 
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U), width = 0.01)+
  theme_bw() +
  facet_wrap(factor(variable,
                    c('aic','bic','ebic'),
                    c('AIC', 'BIC', 'EBIC'))~., 
             scales = 'free',
             nrow = 1)+
  theme(legend.position = 'none') +
  labs(x = '\n spike standard deviation', y = '') +
  scale_y_continuous(expand = c(0.01,0.01), trans='log10') +
  scale_color_manual(values = c('black','grey')) 

ggsave('model_selection_criteria.pdf', width = 9, height = 2.5)

#
ggplot(df_model_criteria_summary[grepl('auc',variable)| variable == 'aic',], aes(x = round(v0,3), y = M, color = type)) + 
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U), width = 0.01)+
  theme_bw() +
  facet_wrap(factor(variable,
                    c('aic','pauc','paucc'),
                    c('AIC', 'pAUC (edge)', 'pAUC (auxiliary variable)'))~.,
             scales = 'free',
             nrow = 1)+
  theme(legend.position = 'none') +
  labs(x = '\n spike standard deviation', y = '') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = c('black','grey')) 

ggsave('model_selection_main.pdf', width = 9, height = 2.5)

########################

### zeta sensitivity ###

########################

dirname <-  '~/simulation_n_200_p_100_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/'
setwd(dirname)

count <- 0
list_GMSS <- list_true_GMSS <- list_GMSS_c <- list_true_GMSS_c <- 
  pt <- index <- it <- vb_it <- vector('list', 9)

for(mn in c(25, 50, 150)){
  for(stddev in c(25, 50, 150)){
    
    count <- count + 1
    list_GMSS[[count]] <- list_true_GMSS[[count]] <- list_GMSS_c[[count]] <- list_true_GMSS_c[[count]] <-
      pt[[count]] <- index[[count]]<- it[[count]] <- vb_it[[count]]<- vector('list', 32)
    
    for(seed in 1:32){
      #
      load(paste0('out_GMSS_VBECM_seed_',seed,'_zeta_mean_',mn,'_sd_',stddev,'.rda'))
      bool_up <- upper.tri(net$A)
      #
      list_GMSS[[count]][[seed]] <- out$estimates$m_delta[bool_up]
      list_true_GMSS[[count]][[seed]] <- net$A[bool_up]
      list_GMSS_c[[count]][[seed]] <- out$estimates$m_gamma
      list_true_GMSS_c[[count]][[seed]] <- beta_true!=0
      #
      tmp <- out$total_pt
      units(tmp) <- 'secs'
      pt[[count]][[seed]] <- tmp
      index[[count]][[seed]] <- out$index
      it[[count]][[seed]] <- sapply(out$full_output, function(x)x$it)
      vb_it[[count]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    }
  }
}

cbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288")

pdf('q50q03_proc_zeta_hyper.pdf', width = 7 ,height = 4)

par(mfrow=c(1,2), pty='s')

# proc edge 
#
for(count in 1:length(list_GMSS)){
  
  plot_roc(list_GMSS[[count]], list_true_GMSS[[count]], col = cbPalette[(count-1)%/% 3 + 1], 
           lty = (count-1)%%3 + 1, add = ifelse(count==1, FALSE, TRUE), fpr_stop = 0.1, 
           main = 'pROC curves \n (edge selection)')
  
}

legend("bottomright",            
       legend = c('5e-3', '1e-2', '3e-2','5e-3', '1e-2', '3e-2'),
       lty = c(rep(1,3), 1:3),
       cex = 0.5,
       col = c(cbPalette[1:3], rep('black',3)), 
       title = 'sparsity (mean & sd)',ncol=2)

# proc auxiliary variable
#
for(count in 1:length(list_GMSS)){
  
  plot_roc(list_GMSS_c[[count]], list_true_GMSS_c[[count]], col = cbPalette[(count-1)%/% 3 + 1], 
           lty = (count-1)%%3 + 1, add = ifelse(count==1, FALSE, TRUE), fpr_stop = 0.1, 
           main = 'pROC curves \n (auxiliary variable selection)')
  
}

legend("bottomright",            
       legend = c('5e-3', '1e-2', '3e-2','5e-3', '1e-2', '3e-2'),
       lty = c(rep(1,3), 1:3),
       cex = 0.5,
       col = c(cbPalette[1:3], rep('black',3)), 
       title = 'sparsity (mean & sd)',ncol=2)

dev.off()

# PPI auxiliary variables
#
m_gamma <- lapply(list_GMSS_c,function(x)do.call('rbind',x))
df_m_gamma <- lapply(m_gamma, function(x)data.table(
  M=apply(x, 2, mean),
  SD=apply(x, 2, sd),
  N=nrow(x),
  ID = seq_len(ncol(x))
))
names(df_m_gamma) <- paste0(rep(c('5e-3', '1e-2', '3e-2'), each = 3), '_',
                            rep(c('5e-3', '1e-2', '3e-2'), times = 3))
df_m_gamma <- rbindlist(df_m_gamma, use.names = T,idcol = 'pars')
df_m_gamma$spar <- as.numeric(gsub('_.*$','',df_m_gamma$pars))
df_m_gamma$sd <- as.numeric(gsub('^.*_','',df_m_gamma$pars))
df_m_gamma$sd <- paste0('sd = ',df_m_gamma$sd)
df_m_gamma$spar <- paste0('mean = ',df_m_gamma$spar)

df_m_gamma[,L:=pmax(M-1.96*SD/sqrt(N), 0)]
df_m_gamma[,U:=pmin(M+1.96*SD/sqrt(N),1)]
tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = beta_true!=0)
df_m_gamma <- merge(df_m_gamma, tmp, by=c('ID'))
df_m_gamma$non_zero <- factor(df_m_gamma$non_zero, c(FALSE,TRUE), c('no','yes'))

ggplot(df_m_gamma, aes(x=ID, y=M, color = non_zero)) + 
  geom_point()+
  geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
  labs(x='\n auxiliary variables', y = "posterior inclusion probability \n",
       color = 'active', title = 'Auxiliary variable PPIs \n')+
  theme_bw() +
  scale_x_discrete()+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  scale_color_grey(start = 0.8, end = 0.2)+
  theme(legend.position = 'bottom',panel.spacing = unit(1, "lines"),
        plot.title =  element_text(hjust = 0.5, face = "bold",size = 20))+
  geom_hline(aes(yintercept = 0.5), linetype = 2) +
  facet_grid(sd~spar) 

ggsave('q50q03_PPI_zeta_hyper.pdf', width = 10, height = 6)




