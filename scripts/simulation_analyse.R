library(navigm)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(xtable)

nrep <- 100
dir_plot <- "~/plots_navigm_out_revised/"

# ###############################
# 
# #######   check files   #######
# 
# ###############################
# 
# cat(" ======== check files ======== \n")
# cat("finish? \n")
# dir_lists <- list.dirs("/rds-d6/user/xx825/hpc-work/navigm_out_revised/",full.names = T)
# dir_lists <- grep("simulation",dir_lists, value = T)
# for (dir in dir_lists) {
#   cat(dir, "\n")
#   cat(length(list.files(dir, pattern = "^out_GM_")), " GM outputs \n")
#   cat(length(list.files(dir, pattern = "^out_GMN_")), " GMN outputs \n")
#   cat(length(list.files(dir, pattern = "^out_GMSS_")), " GMSS outputs \n")
# }
# 
# cat("null graph? \n")
# dir_lists <- list.dirs("/rds-d6/user/xx825/hpc-work/navigm_out_revised/",full.names = T)
# dir_lists <- grep("simulation",dir_lists, value = T)
# for (dir in dir_lists) {
#   cat(dir, "\n")
#   setwd(dir)
#   nedge <- c()
#   for(seed in 1:nrep){
#     load(paste0("data_seed",seed, ".rda"))
#     bool_up <- upper.tri(net$A)
#     nedge <- c(nedge, sum(net$A[bool_up]))
#   }
# 
#   cat(sum(nedge == 0), " graphs have no simulated edges \n")
# }
# cat(" ============================== \n")
# 
# 
###############################

### low number of variables ###

###############################

cat(" =========== low number of variables ============ \n")

cat("Load outputs ... \n")
dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_3_q0_3_zeta_-1.54_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/"
setwd(dirname)

# all replicates
#
list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector("list", nrep)
list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector("list", nrep)

for(seed in 1:nrep){
  
  # gm
  #
  load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  bool_up <- upper.tri(net$A)
  list_GM[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GM[[seed]] <- net$A[bool_up]
  
  # gmn
  #
  load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMN[[seed]] <- net$A[bool_up]
  
  # gmss
  #
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMSS[[seed]] <- net$A[bool_up]
  list_GMSS_c[[seed]] <- out$estimates$m_gamma
  list_true_GMSS_c[[seed]] <- beta_true!=0
}


# compare pROC curves using GM* and GMN
#
cat("Produce pROC curves ... \n")

pdf(file.path(dir_plot, "q3_proc_GM_GMN_GMSS.pdf"), width = 4, height = 4)
par(pty = "s")
plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = "black",
         spread.scale = 1, main = "pROC curves \n (edge selection)",lty = 3)
plot_roc(list_GMN, list_true_GMN, add = T, fpr_stop = 0.1, col = "#818181", spread.scale = 1, lty = 2)
# plot_roc(list_GMSS, list_true_GMSS, add = T, fpr_stop = 0.1, col = "#D4D4D4", spread.scale = 1, lty = 1)
legend("bottomright",
       legend = c("GM*","GMN"),
       lty = c(3,2),
       cex = 0.7,
       col = c("black","#818181"))
dev.off()


# one replicate
#
cat("Check one replicate ... \n")
seed <- 1
load(paste0("data_seed",seed,".rda"))

df_true <- reshape::melt(net$A)
# df_margin_hub <- data.table(margin1 = NA, margin2 = NA, X = 1:nrow(net$A))
# df_margin_hub$margin1 <- V %*% beta_true

#
load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
# set to NULL otherwise https://stackoverflow.com/questions/69666867/constant-warning-message-with-reshapemelt-in-r
rownames(out$estimates$m_delta) <- colnames(out$estimates$m_delta) <- NULL
df <- reshape::melt(out$estimates$m_delta)

#
load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
rownames(out$estimates$m_delta) <- colnames(out$estimates$m_delta) <- NULL
tmp <- reshape::melt(out$estimates$m_delta)

#
beta_est <- matrix(out$estimates$mu_beta, ncol=1)
df_margin_hub <- data.table(margin = NA, X = 1:nrow(net$A))
df_margin_hub$margin <- V %*% beta_est

# gather GM & GMN edge estimation
#
df <- merge(df,tmp, by =c("X1","X2"))
df <- merge(df,df_true, by =c("X1","X2"))
df$GM <- df$GMN <- df$value.plot <- NA
df[df$value.x >= 0.5 & df$value == 1,]$GM <- "TP"
df[df$value.x < 0.5 & df$value == 0,]$GM <- "TN"
df[df$value.y >= 0.5 & df$value == 1,]$GMN <- "TP"
df[df$value.y < 0.5 & df$value == 0,]$GMN <- "TN"
df[is.na(df$GM),]$GM <- "F"
df[is.na(df$GMN),]$GMN <- "F"
df[df$X1 > df$X2,]$value.plot <- df[df$X1 > df$X2,]$GMN
df[df$X1 < df$X2,]$value.plot <- df[df$X1 < df$X2,]$GM
df[df$X1 == df$X2,]$value.plot <- NA

p <- ggplot(df, aes(x=X1, y=X2, fill=value.plot)) +
  geom_tile() +
  scale_fill_grey(start = 0.8, end = 0.2, na.value = "white",na.translate = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(fill = "", x= "", y = "") +
  annotate("text", x = 25, y = 75, label = "GM*", fontface = "bold") +
  annotate("text", x = 75, y = 25, label = "GMN", fontface = "bold") +
  theme(legend.position = "top",
        legend.key.height= unit(0.3, "cm"),
        legend.key.width= unit(0.3, "cm"))

p_legend <- ggpubr::get_legend(p)
p <- p + theme(legend.position = "none")
px <- ggplot(df_margin_hub, aes(X, 1, fill = margin)) +
  geom_tile() +
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_gradient(low = "white", high = "black")+
  labs(fill = "hub propensity") +
  theme(legend.position="right",
        legend.direction = "horizontal",
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(color="white"),
        axis.ticks.x=element_line(color="white"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8))+
  guides(fill = guide_colourbar(barwidth = 5,  barheight = 0.3, title.position = "top"))

px_legend <- ggpubr::get_legend(px)
px <- px + theme(legend.position = "none")

#
df <- setDT(df)
df_margin_sum_tp <- df[X1 != X2, list(sum_GM = sum(GM == "TP"),
                                      sum_GMN = sum(GMN == "TP")), by = "X1"]
df_margin_sum_tn <- df[X1 != X2, list(sum_GM = sum(GM == "TN"),
                                      sum_GMN = sum(GMN == "TN")), by = "X1"]
df_margin_sum_tp  <- melt(df_margin_sum_tp, id.vars = "X1")
df_margin_sum_tn  <- melt(df_margin_sum_tn, id.vars = "X1")

py_tp <-
  ggplot(df_margin_sum_tp, aes(X1, value, fill = factor(variable,c("sum_GM","sum_GMN"), c("GM*", "GMN"))))+
  coord_flip() +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_grey() +
  theme_classic() +
  labs(y = "# TP", fill ="") +
  theme(legend.position="top",
        legend.direction = "vertical",
        legend.box = "vertical",
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.height= unit(0.3, "cm"),
        legend.key.width= unit(0.3, "cm"))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

py_tp_legend <- ggpubr::get_legend(py_tp)
py_tp <- py_tp + theme(legend.position = "none")

p + px + py_tp + p_legend + px_legend + py_tp_legend +
  plot_layout(design = "46
                        13
                        25",
              widths = c(0.8, 0.2),
              heights = c(0.05, 0.85, 0.05))

ggsave(file.path(dir_plot, paste0("q3_heatmap_GM_GMN_seed",seed,".pdf")),
       width = 6, height = 5.5)

#################################################

### high number of variables  + null scenario ###

#################################################

cat(" =========== high number of variables  + null scenario ============ \n")

for(label in c("q_50_q0_0", "q_50_q0_3")){
  
  if (label == "q_50_q0_3"){
    cat(" =========== high number of variables ============ \n")
    dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
  }else if(label == "q_50_q0_0"){
    cat(" =========== null scenario ============ \n")
    dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0"
  }
  
  setwd(dirname)
  
  list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector("list",nrep)
  list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector("list",nrep)
  
  list_GM_omega <- list_GMN_omega <- list_GMSS_omega <- vector("list",nrep)
  list_true_GM_omega <- list_true_GMN_omega <- list_true_GMSS_omega <- vector("list",nrep)
  
  cat("Load outputs ... \n")
  
  for(seed in 1:nrep){
    
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    bool_up <- upper.tri(net$A)
    list_GM[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GM[[seed]] <- net$A[bool_up]
    list_GM_omega[[seed]] <- out$estimates$Omega[bool_up]
    list_true_GM_omega[[seed]] <- net$Omega[bool_up]
    
    # gmn
    #
    load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMN[[seed]] <- net$A[bool_up]
    list_GMN_omega[[seed]] <- out$estimates$Omega[bool_up]
    list_true_GMN_omega[[seed]] <- net$Omega[bool_up]
    
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
    list_GMSS_omega[[seed]] <- out$estimates$Omega[bool_up]
    list_true_GMSS_omega[[seed]] <- net$Omega[bool_up]
    
  }
  
  cat("Produce pROC curves ... \n")
  
  pdf(file.path(dir_plot, paste0(gsub("_","",label),"_proc_GM_GMN_GMSS.pdf")),width = 4, height = 4)
  par(pty = "s")
  plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1, main = "pROC curves \n (edge selection)",lty = 3)
  plot_roc(list_GMN, list_true_GMN, add = T, fpr_stop = 0.1, col = "#818181", spread.scale = 1, lty = 2)
  plot_roc(list_GMSS, list_true_GMSS, add = T, fpr_stop = 0.1, col = "#D4D4D4", spread.scale = 1, lty = 1)
  legend("bottomright",
         legend = c("GM*","GMN","GMSS"),
         lty = 3:1,
         cex = 0.7,
         col = c("black","#818181","#D4D4D4"))
  dev.off()
  
  
  cat("Plot PPIs ... \n")
  
  tmp <- do.call("rbind",list_GMSS_c)
  df_m_gamma <- data.table(M = apply(tmp, 2, mean),
                           SD = apply(tmp, 2, sd),
                           N = nrow(tmp),
                           ID = seq_len(ncol(tmp)))
  df_m_gamma[,L:=pmax(M-SD/sqrt(N), 0)]
  df_m_gamma[,U:=pmin(M+SD/sqrt(N), 1)]
  tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = list_true_GMSS_c[[1]]!=0)
  df_m_gamma <- merge(df_m_gamma, tmp, by=c("ID"))
  ggplot(df_m_gamma, aes(x=ID, y=M, color=non_zero)) +
    geom_point()+
    geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
    labs(x="\n  auxiliary variables", y = "posterior inclusion probability \n",
         color = "active", title = "Auxiliary variable PPIs\n")+
    theme_bw() +
    scale_x_discrete()+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))+
    scale_color_manual(breaks = c(FALSE,TRUE),
                       labels = c("no","yes"),
                       values=c("#BABABA","black"))+
    theme(legend.position = "bottom",panel.spacing = unit(1, "lines"))+
    geom_hline(aes(yintercept = 0.5), linetype = 2) +
    theme(plot.title =  element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(dir_plot,paste0(gsub("_","",label),"_PPI_GMSS_c.pdf")), width = 4, height = 3.5)
  
  cat(" ======================================== \n")
}


##########################

### reference scenario ###

##########################

cat(" =========== reference scenario ============ \n")
#
cat("Edges inferred by GM* or GMSS: \n")
num_in_gm_not_in_gmss <- num_in_gmss_not_in_gm <- c()
true_in_gm_not_in_gmss <- true_in_gmss_not_in_gm <- list()
dfppi_in_gmss_not_in_gm <- dfppi_in_gm_not_in_gmss <- list()

for (seed in 1:nrep){
  #
  num_in_gm_not_in_gmss  <- c(num_in_gm_not_in_gmss, sum(list_GM[[seed]] > 0.5 & list_GMSS[[seed]]< 0.5))
  num_in_gmss_not_in_gm <- c(num_in_gmss_not_in_gm, sum(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]>0.5))
  true_in_gm_not_in_gmss[[seed]] <- list_true_GM[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]]< 0.5)]
  true_in_gmss_not_in_gm[[seed]] <- list_true_GM[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]> 0.5)]
  #
  if(num_in_gmss_not_in_gm[seed]!=0){
    tmp <- data.table(`PPI_GM*` = list_GM[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]>0.5)],
                      PPI_GMSS = list_GMSS[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]>0.5)],
                      EDGE_TRUE = list_true_GM[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]> 0.5)],
                      `Omega_GM*` = list_GM_omega[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]> 0.5)],
                      Omega_GMSS = list_GMSS_omega[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]> 0.5)],
                      Omega_TRUE = list_true_GM_omega[[seed]][which(list_GM[[seed]] < 0.5 & list_GMSS[[seed]]> 0.5)],                      
                      ID = 1:num_in_gmss_not_in_gm[seed],
                      SEED = seed)
    dfppi_in_gmss_not_in_gm[[seed]] <- tmp
  }
  #
  if(num_in_gm_not_in_gmss[seed]!=0){
    tmp <- data.table(`PPI_GM*` = list_GM[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],
                      PPI_GMSS = list_GMSS[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],
                      EDGE_TRUE = list_true_GM[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],
                      `Omega_GM*` = list_GM_omega[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],
                      Omega_GMSS = list_GMSS_omega[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],
                      Omega_TRUE = list_true_GM_omega[[seed]][which(list_GM[[seed]] > 0.5 & list_GMSS[[seed]] < 0.5)],                      
                      ID = 1:num_in_gm_not_in_gmss[seed],
                      SEED = seed)
    dfppi_in_gm_not_in_gmss[[seed]] <- tmp
  }
}

#
dfppi_in_gmss_not_in_gm <- rbindlist(dfppi_in_gmss_not_in_gm)
dfppi_in_gm_not_in_gmss <- rbindlist(dfppi_in_gm_not_in_gmss)
tmp <- rbind(dfppi_in_gmss_not_in_gm[EDGE_TRUE == 1,],dfppi_in_gm_not_in_gmss[EDGE_TRUE==1, ])
tmp[,`ABS_ERR_GM*`:=abs(`Omega_GM*`-Omega_TRUE)]
tmp[,ABS_ERR_GMSS:=abs(Omega_GMSS-Omega_TRUE)]
tmp <- tmp[,list(`ABS_ERR_GM*` = mean(`ABS_ERR_GM*`),
                 ABS_ERR_GMSS = mean(ABS_ERR_GMSS)), by = "SEED"]

cat("MAE (GM*): \n")
print(paste0(round(mean(tmp$`ABS_ERR_GM*`),2), "(",round(sd(tmp$`ABS_ERR_GM*`)/sqrt(nrep),2),")"))
cat("MAE (GMSS): \n")
paste0(round(mean(tmp$`ABS_ERR_GMSS`),2), "(",round(sd(tmp$`ABS_ERR_GMSS`)/sqrt(nrep),2), ")")

#
cat("In GMSS not in GM* PPIs (GM* and GMSS resp): \n")
print(range(dfppi_in_gmss_not_in_gm$`PPI_GM*`))
print(range(dfppi_in_gmss_not_in_gm$`PPI_GMSS`))
cat("In GM* not in GMSS PPIs (GM* and GMSS resp): \n")
print(range(dfppi_in_gm_not_in_gmss$`PPI_GM*`))
print(range(dfppi_in_gm_not_in_gmss$`PPI_GMSS`))

#
tmp <- rbind(data.table(prop = sapply(true_in_gm_not_in_gmss, function(x)sum(x)/length(x)),
                        n = sapply(true_in_gm_not_in_gmss,length),
                        id = 1:length(true_in_gm_not_in_gmss),
                        type = "in GM* not in GMSS"),
             data.table(prop = sapply(true_in_gmss_not_in_gm, function(x)sum(x)/length(x)),
                        n = sapply(true_in_gmss_not_in_gm,length),
                        id = 1:length(true_in_gmss_not_in_gm),
                        type = "in GMSS not in GM*"))

#
ggplot(tmp,aes(id,prop,color = type, size = n))+
  geom_point(position = position_dodge(width=0.05), alpha = 0.6) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal") +
  guides(size = guide_legend(nrow = 1,override.aes = list(shape=1))) +
  labs(x = "replicates", y = "proportions of true positives among edges \n inferred in either GM* or GMSS", color = "edges reported: ", size = "number of edges") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(limits = c(-1,101), breaks = seq(0,100,by =5)) +
  scale_color_grey() 
ggsave(file.path(dir_plot,"scatter_unique_edges.pdf"), width = 10, height=4)

#
df <- data.table(omega = list_GMSS_omega[[seed]],
                 delta = list_GMSS[[seed]] > 0.5)
df$delta <- factor(df$delta, c(T,F), c("edges inferred as present","edges inferred as absent"))
ggplot(df, aes(omega)) +
  geom_histogram() +
  facet_wrap(delta~., scales = "free", nrow = 1) +
  labs(x="\n precision matrix entry estimates",fill = "") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(dir_plot, paste0("q50q03_omega_hist_seed_",seed,".pdf")),
       width = 8, height = 3)

#####################

### null scenario ###

#####################

cat(" =========== null scenario ============ \n")
# estimated effects plot
#

dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0"

setwd(dirname)
m1_beta_GMN <- m1_beta_GMSS <-  m1_zeta_GMN <- m1_zeta_GMSS <- list()
pauc_GMN <- pauc_GMSS <- c()

for(seed in 1:nrep){
  
  load(paste0("data_seed",seed,".rda"))
  bool_up <- upper.tri(net$A)
  
  # gmn
  #
  load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  m1_beta_GMN[[seed]] <- out$estimates$mu_beta * out$estimates$m_gamma
  m1_zeta_GMN[[seed]] <- out$estimates$mu_zeta
  pauc_GMN <- c(pauc_GMN, compute_pauc(out$estimates$m_delta[bool_up], net$A[bool_up], fpr_stop = 0.1, standardise = T))
  
  # gmss
  #
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  m1_beta_GMSS[[seed]] <- out$estimates$mu_beta * out$estimates$m_gamma
  m1_zeta_GMSS[[seed]] <- out$estimates$mu_zeta
  pauc_GMSS <- c(pauc_GMSS, compute_pauc(out$estimates$m_delta[bool_up], net$A[bool_up], fpr_stop = 0.1, standardise = T))
}

print("pAUC (GMN) ..")
print(mean(pauc_GMN))
print(sd(pauc_GMN)/sqrt(nrep))
print("pAUC (GMSS) ..")
print(mean(pauc_GMSS))
print(sd(pauc_GMSS)/sqrt(nrep))


print("zeta (GMN) ..")
print(quantile(do.call("c",m1_zeta_GMN)))
print("zeta (GMSS) ..")
print(quantile(do.call("c",m1_zeta_GMSS)))

pdf(file.path(dir_plot,"q50q00_zeta_gmn_gmss.pdf"),width = 6, height = 4)
plot(do.call("c",m1_zeta_GMN),do.call("c",m1_zeta_GMSS),xlab = "zeta (GMN)", ylab = "zeta (GMSS)")
abline(a=0,b=1, lty=2)
dev.off()

cat("beta .. \n")
m1_beta_GMN <- do.call("rbind", m1_beta_GMN)
m1_beta_GMSS <- do.call("rbind", m1_beta_GMSS)

df <- rbind(
  data.table(M = apply(m1_beta_GMN,2,mean),
             SD = apply(m1_beta_GMN,2,sd),
             N = apply(m1_beta_GMN,2,length),
             ID = seq_len(ncol(m1_beta_GMN)),
             method = "GMN-VBECM"),
  data.table(M = apply(m1_beta_GMSS,2,mean),
             SD = apply(m1_beta_GMSS,2,sd),
             N = apply(m1_beta_GMSS,2,length),
             ID = seq_len(ncol(m1_beta_GMSS)),
             method = "GMSS-VBECM")
)

df[,L:=M-SD/sqrt(N)]
df[,U:=M+SD/sqrt(N)]

# dfs <- dcast(df, ID ~ method,value.var = c("M","L","U"))
# dfs$active <- beta_true!=0
# tmp <- range(c(df$M, df$L, df$U))
# 
# ggplot(dfs, aes(`M_GMSS-VBECM` , `M_GMN-VBECM`))+
#   geom_errorbar(mapping = aes(ymin=`L_GMN-VBECM`, ymax=`U_GMN-VBECM`), color = "grey") +
#   geom_errorbarh(mapping = aes(xmin=`L_GMSS-VBECM`,xmax=`U_GMSS-VBECM`), color = "grey") +
#   geom_point() +
#   theme_classic() +
#   theme(legend.position = "bottom") +
#   theme(plot.title =  element_text(hjust = 0.5, face = "bold"))+
#   labs(x="estimated effects GMSS", y="estimated effects GMN",
#        title = "Auxiliary variable effects\n") +
#   coord_fixed()+
#   scale_x_continuous(limits=tmp) +
#   scale_y_continuous(limits=tmp) +
#   geom_hline(yintercept = 0, linetype=2)+
#   geom_vline(xintercept = 0, linetype=2)
# ggsave(file.path(dir_plot,"q50q00_scatter_effects_GMSS_GMN_c.pdf"),
#        width = 3, height = 3)


ggplot(df, aes(ID , M, color = method))+
  geom_errorbar(mapping = aes(ymin = L, ymax=U, color = method)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="auxiliary variables", y="estimated effects", color = "")+
  scale_color_grey()
ggsave(file.path(dir_plot,"q50q00_effects_GMSS_GMN_c.pdf"),
       width = 6, height = 4.5)


tmp <- rbind(data.table(value = m1_beta_GMN[1,],
                        id = 1: length(m1_beta_GMN[1,]),
                        method = "GMN"),
             data.table(value = m1_beta_GMSS[1,],
                        id = 1: length(m1_beta_GMSS[1,]),
                        method = "GMSS"))

ggplot(tmp, aes(value, fill = factor(method, c("GMSS","GMN")))) +
  geom_histogram(alpha = 0.5, position = "identity") +
  theme_bw() +
  theme(legend.position = "bottom") + 
  scale_x_continuous(limits = c(-0.55,0.55), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,51), expand = c(0,0)) +
  labs(x="estimated effects", fill = "")+
  scale_fill_grey(start = 0.8, end = 0.2) 
ggsave(file.path(dir_plot,"q50q00_effects_GMSS_GMN_c_data1.pdf"),
       width = 4, height = 4)


#########################

### inference methods ###

#########################

cat(" =========== inference methods ============ \n")

# GM
#
cat("GM: \n")
dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_0_zeta_-1.85_noise_0.1_beta0_0.5_codata_FALSE_hub_FALSE_corr_0/"
setwd(dirname)

list_GM <- list_true_GM <- list_GM_ECM <- list_true_GM_ECM <- list()
pt_vbem <- index_vbem  <- it_vbem  <- vb_it_vbem <- list()
pt_em <- index_em  <- it_em  <- vb_it_em <- list()

cat("Load outputs ...\n")
for(seed in 1:nrep){
  
  # vbem
  #
  load(paste0("out_GM_v1_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01.rda"))
  bool_up <- upper.tri(out$estimates$m_delta)
  list_GM[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GM[[seed]] <- net$A[bool_up]
  
  #
  tmp <- out$total_pt
  units(tmp) <- "secs"
  pt_vbem[[seed]] <- tmp
  index_vbem[[seed]] <- out$index
  it_vbem[[seed]] <- sapply(out$full_output, function(x)x$it)
  vb_it_vbem[[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it)) # max iterations in vb
  
  # em
  #
  load(paste0("out_GM_v1_ECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01.rda"))
  list_GM_ECM[[seed]] <- out$estimates$P1[bool_up]
  list_true_GM_ECM[[seed]] <- net$A[bool_up]
  
  #
  tmp <- out$total_pt
  units(tmp) <- "secs"
  pt_em[[seed]] <- tmp
  index_em[[seed]] <- out$index
  it_em[[seed]] <- sapply(out$full_output, function(x)x$it)
  
}

cat("Plot pROC curves ...\n")

pdf(file.path(dir_plot, "proc_GM_em_vbem.pdf"), width = 4, height = 4)
par(mfrow = c(1,1), pty = "s")
plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = "black", spread.scale = 1, main="pROC curves \n (edge selection)",lty=2)
plot_roc(list_GM_ECM, list_true_GM_ECM, add = T, fpr_stop = 0.1,  spread.scale = 1, col = "grey")
legend("bottomright",
       legend = c("GM-VBECM","GM-ECM"),
       lty = c(2,1),
       cex = 0.7,
       col = c("black","grey"))
dev.off()

# runtime
#
cat("Runtime ...\n")

print("runtime (ecm)")
print(round(mean(unlist(pt_em)),2))
print(round(sd(unlist(pt_em))/length(unlist(pt_em)),2))

print("runtime (vbecm)")
print(round(mean(unlist(pt_vbem)),2))
print(round(sd(unlist(pt_vbem))/length(unlist(pt_vbem)),2))

# GMSS
#
cat("GMSS: \n")
dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/"
setwd(dirname)

list_GMSS <- list_GMSS_c  <- list_true_GMSS <- list_true_GMSS_c <- vector("list",nrep)
list_GMSS_ECM <- list_GMSS_ECM_c <- list_true_GMSS_ECM <- list_true_GMSS_ECM_c <- vector("list",nrep)
pt_em <- index_em <- it_em <- vb_it_em <- list()
pt_vbem <- index_vbem <- it_vbem <- vb_it_vbem <- list()

for(seed in 1:nrep){
  
  # vbem
  #
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01.rda"))
  bool_up <- upper.tri(out$estimates$m_delta)
  list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
  list_true_GMSS[[seed]] <- net$A[bool_up]
  list_GMSS_c[[seed]] <- out$estimates$m_gamma
  list_true_GMSS_c[[seed]] <- beta_true!=0
  
  #
  tmp <- out$total_pt
  units(tmp) <- "secs"
  pt_vbem[[seed]] <- tmp
  index_vbem[[seed]] <- out$index
  it_vbem[[seed]] <- sapply(out$full_output, function(x)x$it)
  vb_it_vbem[[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
  
  # em
  #
  load(paste0("out_GMSS_ECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01.rda"))
  list_GMSS_ECM[[seed]] <- out$estimates$P1[bool_up]
  list_true_GMSS_ECM[[seed]] <- net$A[bool_up]
  list_GMSS_ECM_c[[seed]] <- out$estimates$P2
  list_true_GMSS_ECM_c[[seed]] <- beta_true!=0
  
  #
  tmp <- out$total_pt
  units(tmp) <- "secs"
  pt_em[[seed]] <- tmp
  index_em[[seed]] <- out$index
  it_em[[seed]] <- sapply(out$full_output, function(x)sapply(x,function(y)y$it))
  
}

cat("plot pROC curves .. \n")
pdf(file.path(dir_plot, "proc_GMSS_em_vbem.pdf"), width = 4, height = 4)
par(mfrow = c(1,1), pty = "s")
# edges
#
plot_roc(list_GMSS, list_true_GMSS, add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1,
         main = "pROC curves \n (edge selection)", lty=2, ylab = "")
plot_roc(list_GMSS_ECM, list_true_GMSS_ECM, add = T, fpr_stop = 0.1,  spread.scale = 1, col = "grey")
legend("bottomright",
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2,1),
       cex = 0.7,
       col = c("black","grey"))
dev.off()

# auxiliary variables
#
pdf(file.path(dir_plot, "procc_GM_em_vbem.pdf"), width = 4, height = 4)
par(mfrow = c(1,1), pty = "s")
plot_roc(list_GMSS_c, list_true_GMSS_c, add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1,
         main = "pROC curves \n (auxiliary variable selection)", ylab = "",lty = 2)
plot_roc(list_GMSS_ECM_c, list_true_GMSS_ECM_c, add = T, fpr_stop = 0.1,  spread.scale = 1, col = "grey")
legend("bottomright",
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2, 1),
       cex = 0.7,
       col = c("black","grey"))
dev.off()

# runtime
#
cat("Runtime ..\n")
print("runtime (ecm)")
print(round(mean(unlist(pt_em)), 2)/60)
print(round(sd(unlist(pt_em))/length(unlist(pt_em)), 2))
print("runtime (vbecm)")
print(round(mean(unlist(pt_vbem)), 2))
print(round(sd(unlist(pt_vbem))/length(unlist(pt_vbem)), 2))

############################

### sensitivity analysis ###

############################

cat(" =========== sensitivity analysis ============ \n")

vec_dirnames <- c(
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_50_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_100_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_100_p_50_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_20_q0_3_zeta_-1.53_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_100_q0_3_zeta_-1.54_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_1_zeta_-0.79_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_5_zeta_-2.16_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.55_noise_0.2_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.56_noise_0.3_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.8_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0",
  "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.21_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
)


fpr_stop <- 0.1
cat("fpr_stop = ",fpr_stop, "...\n")
pauc_GMSS <- pauc_GMSS_c <-  pauc_GMN <-  pauc_GM <- list()
pt_GM <- index_GM <- it_GM <- vb_it_GM <- list()
pt_GMN <- index_GMN <- it_GMN <- vb_it_GMN  <- list()
pt_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS  <- list()

for (dirname in vec_dirnames){
  
  cat(dirname, ".. \n")
  setwd(dirname)
  
  list_GMSS <- list_true_GMSS <- list_GMN <- list_true_GMN <- list_GM <- list_true_GM <- list()
  list_GMSS_c <- list_true_GMSS_c <- list()
  
  pt_GM[[dirname]] <- index_GM[[dirname]]<- it_GM[[dirname]] <- vb_it_GM[[dirname]] <- list()
  pt_GMN[[dirname]] <- index_GMN[[dirname]]<- it_GMN[[dirname]] <- vb_it_GMN[[dirname]] <- list()
  pt_GMSS[[dirname]] <- index_GMSS[[dirname]]<- it_GMSS[[dirname]] <- vb_it_GMSS[[dirname]] <- list()
  
  for(seed in 1:nrep){
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    bool_up <- upper.tri(out$estimates$m_delta)
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
    
    tmp <- out$total_pt
    units(tmp) <- "secs"
    pt_GMSS[[dirname]][[seed]] <- tmp
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMN[[seed]] <- net$A[bool_up]
    
    tmp <- out$total_pt
    units(tmp) <- "secs"
    pt_GMN[[dirname]][[seed]] <- tmp
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    list_GM[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GM[[seed]] <- net$A[bool_up]
    
    tmp <- out$total_pt
    units(tmp) <- "secs"
    pt_GM[[dirname]][[seed]] <- tmp
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
  id <- which(nedge == 0)
  if(length(id) != 0){
    pauc_GMSS[[dirname]] <- compute_pauc(list_GMSS[-id], list_true_GMSS[-id], fpr_stop = fpr_stop, standardise = T)
    pauc_GMN[[dirname]] <- compute_pauc(list_GMN[-id], list_true_GMN[-id], fpr_stop = fpr_stop, standardise = T)
    pauc_GM[[dirname]] <- compute_pauc(list_GM[-id], list_true_GM[-id], fpr_stop = fpr_stop, standardise = T)
  }else{
    pauc_GMSS[[dirname]] <- compute_pauc(list_GMSS, list_true_GMSS, fpr_stop = fpr_stop, standardise = T)
    pauc_GMN[[dirname]] <- compute_pauc(list_GMN, list_true_GMN, fpr_stop = fpr_stop, standardise = T)
    pauc_GM[[dirname]] <- compute_pauc(list_GM, list_true_GM, fpr_stop = fpr_stop, standardise = T)
  }
  
  pauc_GMSS_c[[dirname]] <- compute_pauc(list_GMSS_c, list_true_GMSS_c, fpr_stop = fpr_stop, standardise =T)
}

# check index and iterations
print("Check indices and iterations ..")
print("GM:")
print(sapply(lapply(index_GM, unlist), quantile))
print(sapply(lapply(index_GM, unlist), function(x)sum(x==1)))
print(sapply(lapply(index_GM, unlist), function(x)sum(x==16)))
print(any(sapply(lapply(it_GM, unlist), max) == 1e5))
print(any(sapply(lapply(vb_it_GM, unlist), max)== 1e5))

print("GMN:")
print(sapply(lapply(index_GMN, unlist), quantile))
print(sapply(lapply(index_GMN, unlist), function(x)sum(x==1)))
print(sapply(lapply(index_GMN, unlist), function(x)sum(x==16)))
print(any(sapply(lapply(it_GMN, unlist), max) == 1e5))
print(any(sapply(lapply(vb_it_GMN, unlist), max)== 1e5))

print("GMSS:")
print(sapply(lapply(index_GMSS, unlist), quantile))
print(sapply(lapply(index_GMSS, unlist), function(x)sum(x==1)))
print(sapply(lapply(index_GMSS, unlist), function(x)sum(x==16)))
print(any(sapply(lapply(it_GMSS, unlist), max) == 1e5))
print(any(sapply(lapply(vb_it_GMSS, unlist), max)== 1e5))

#
# create table for pauc
#
print("Create pAUC tables ..")
pauc_GMSS <- do.call("rbind",  pauc_GMSS)
pauc_GMSS <- apply(pauc_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GMSS_c <- do.call("rbind",  pauc_GMSS_c)
pauc_GMSS_c <- apply(pauc_GMSS_c, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GMN <- do.call("rbind",  pauc_GMN)
pauc_GMN <- apply(pauc_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
pauc_GM <- do.call("rbind",  pauc_GM)
pauc_GM <- apply(pauc_GM, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(pauc_GM),t(pauc_GMN),t(pauc_GMSS), t(pauc_GMSS_c)), keep.rownames = T)
df$variable <- rep(c("pauc_GM","pauc_GMN","pauc_GMSS","pauc_GMSS_c"), each = length(unique(df$rn)))
setnames(df, c(paste0("V",1:3)), c("M","SD","N"))
df[,label:=paste0(round(M,2), " (", round(SD/sqrt(N), 2)  ,") ")]
tmp <- dcast(df, rn~variable, value.var = "label")
tmp[,rn := factor(rn, vec_dirnames)]

setkey(tmp, rn)
print(xtable(tmp[,c("rn","pauc_GM","pauc_GMN","pauc_GMSS","pauc_GMSS_c")]), include.rownames=FALSE)

print("Create runtime tables ..")
#
# create table for runtime
#
pt_GMSS_v <- do.call("rbind",  lapply(pt_GMSS,"unlist"))
pt_GMSS_v <- apply(pt_GMSS_v, 1, function(x){c(mean(x), sd(x), length(x))})
pt_GMN_v <- do.call("rbind",  lapply(pt_GMN,"unlist"))
pt_GMN_v <- apply(pt_GMN_v, 1, function(x){c(mean(x), sd(x), length(x))})
pt_GM_v <- do.call("rbind",  lapply(pt_GM,"unlist"))
pt_GM_v <- apply(pt_GM_v, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(pt_GM_v),t(pt_GMN_v),t(pt_GMSS_v)), keep.rownames = T)
df$variable <- rep(c("time_GM","time_GMN","time_GMSS"), each = length(unique(df$rn)))
setnames(df, c(paste0("V",1:3)), c("M","SD","N"))
df[,label:=paste0(round(M,2), " (", round(SD/sqrt(N), 2)  ,") ")]
tmp <- dcast(df, rn~variable, value.var = "label")
tmp[,rn := factor(rn, vec_dirnames)]

setkey(tmp, rn)
print(xtable(tmp[,c("rn","time_GM","time_GMN","time_GMSS")]), include.rownames = F)


#######################

### runtime profile ###

#######################

cat(" =========== runtime analysis ============ \n")

# sample size

cat("Effect of sample sizes: \n")

runtime_GM <- index_GM <- it_GM <- vb_it_GM <- list()
runtime_GMN <- index_GMN <- it_GMN <- vb_it_GMN <- list()
runtime_GMSS <-  index_GMSS <- it_GMSS <- vb_it_GMSS <- list()

cat("Load outputs .. \n")

# for(n in c(20, 50, 100, 200, 300, 400, 500)){
for(n in c(20, 50, 100, 200, 300)){
  
  cat("n = ", n, " \n")
  dirname <-  paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_",n,"_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  
  setwd(dirname)
  
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  for(seed in 1:nrep){
    
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]] <- c(index_GM[[dirname]] ,out$index)
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
    index_GMN[[dirname]][[seed]] <- out$index
    it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

cat("Plot runtime ..\n")
runtime_GM <- do.call("rbind", runtime_GM)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- do.call("rbind", runtime_GMN)
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMSS <- do.call("rbind", runtime_GMSS)
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})


df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c("GM","GMN","GMSS"), each = length(unique(df$rn)))
setnames(df, c(paste0("V",1:3)), c("M","SD","N"))
df[,L:=M- SD/sqrt(N)]
df[,U:=M+ SD/sqrt(N)]
df[,N := as.numeric(str_extract(rn, "(?<=n_)\\d+"))]
# save(df, file = file.path(dir_plot, "runtime_vs_n.rda"))


df[variable=="GM", variable:="GM*"]
ggplot(df, aes(x=N, y=M, group=variable, color = variable, lty=variable)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width= 10) +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x= "sample size (N)", y="second",color ="",lty="")  +
  scale_color_grey() 

ggsave(file.path(dir_plot, "runtime_vs_n.pdf"),
       width = 4, height = 3)

ggplot(df, aes(x=N, y=M, group=variable, color = variable, lty=variable)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width= 10) +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x= "sample size (N)", y="second (log-scale)",color ="",lty="")  +
  scale_y_continuous(trans="log10") +
  scale_color_grey() 

ggsave(file.path(dir_plot, "log_runtime_vs_n.pdf"),
       width = 4, height = 3)


# graph size

cat("Effect of graph sizes: \n")
runtime_GM <- index_GM <- it_GM <- vb_it_GM <- list()
runtime_GMN <- index_GMN <- it_GMN <- vb_it_GMN <- list()
runtime_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS <- list()

for(p in c(50, seq(100, 600, by = 100), 800)){
  cat("p = ", p, " \n")
  if(p == 50){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_",p,"_q_50_q0_3_zeta_-1.5_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }else if(p == 100 | p == 300){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_",p,"_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }else if(p == 200 | p == 600){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_",p,"_q_50_q0_3_zeta_-1.54_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }else if(p == 400 | p == 800 | p == 1000){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_",p,"_q_50_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }else if(p == 500){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_",p,"_q_50_q0_3_zeta_-1.53_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }
  
  setwd(dirname)
  
  #
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  cat("Load outputs ..\n")
  #
  for(seed in 1:nrep){
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    if(file.exists(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))){
      load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
      tmp <- out$total_pt
      units(tmp) <- "secs"
      runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
      index_GMN[[dirname]][[seed]] <- out$index
      it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
      vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    }else{
      runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], NA)
      index_GMN[[dirname]][[seed]] <- NA
      it_GMN[[dirname]][[seed]] <- NA
      vb_it_GMN[[dirname]][[seed]] <- NA
    }
    
    
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

cat("Plot runtime ..\n")
runtime_GM <- do.call("rbind", runtime_GM)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- do.call("rbind", runtime_GMN)
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x,na.rm = T), sd(x, na.rm = T), sum(!is.na(x)))})
runtime_GMSS <- do.call("rbind", runtime_GMSS)
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})

df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c("GM","GMN","GMSS"), each = length(unique(df$rn)))
setnames(df, c(paste0("V",1:3)), c("M","SD","N"))
df[,L:=M-SD/sqrt(N)]
df[,U:=M+SD/sqrt(N)]
df[,P := as.numeric(str_extract(rn, "(?<=p_)\\d+"))]
df <- df[df$N > ceiling(0.9 * nrep),]
#
# print(df[df$variable == "GMSS",]$M/60/60)
# save(df, file = file.path(dir_plot, "runtime_vs_p.rda"))


df[variable=="GM",variable:="GM*"]
ggplot(df, aes(x=P, y=M, group=variable, color = variable, lty=variable)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width = 10) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x= "graph size (P)", y="second ",color ="",lty="")  +
  scale_color_grey()
ggsave(file.path(dir_plot, "runtime_vs_p.pdf"),
       width = 4, height = 3)

ggplot(df, aes(x=P, y=M, group=variable, color = variable, lty=variable)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width = 10) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x= "graph size (P)", y="second (log-scale)",color ="",lty="")  +
  scale_y_continuous(trans="log10")+
  scale_color_grey()
ggsave(file.path(dir_plot, "log_runtime_vs_p.pdf"),
       width = 4, height = 3)


# breakdown
#
dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0/"
setwd(dirname)
pt <-list()

for(seed in 1:nrep){
  #
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
  pt[[seed]] <-
    sapply(out$full_output, function(x){
      tmp <- x$pt
      units(tmp) <- "secs"
      return(tmp)
    })
}

#
pt <- do.call("cbind",pt)
pt <- apply(pt, 1, function(x)c(mean(x),sd(x),length(x)))
pt <- data.table(t(pt))
pt$v0 <- seq(1e-2, 1, length.out = 16)
colnames(pt)[1:3] <- c("M","SD","N")
pt[,L:= M-SD/sqrt(N)]
pt[,U:= M+SD/sqrt(N)]

ggplot(pt, aes(v0, M)) +
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U),width = 0.02) +
  theme_bw() +
  scale_y_continuous(trans = "log10") +
  labs(x="\n spike standard deviation", y = "second (log-scale)") +
  theme(
    plot.margin = margin(5.5,5.5,25.5,5.5)
  )
ggsave(file.path(dir_plot, "q50q03_runtime_GMSS_reference_gridstart_0.01_init0.pdf"),
       width = 4, height = 3)


# number of auxiliary variables
#
cat("Effects of auxiliary variables: \n ")

runtime_GM <-   index_GM <- it_GM <- vb_it_GM <- list()
runtime_GMN <- index_GMN <- it_GMN <- vb_it_GMN <- list()
runtime_GMSS <- index_GMSS <- it_GMSS <- vb_it_GMSS <- list()

cat("Load outputs .. \n")
for(q in c(20, 50, 100, 200, 300, 400, 600, 800, 1000)){
  cat("q = ", q, "\n")
  
  #
  if(q == 20){
    dirname <-"/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_20_q0_3_zeta_-1.53_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
  }else if(q == 50){
    dirname <-"/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
  }else if(q >=100){
    dirname <- paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_",q,"_q0_3_zeta_-1.54_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
  }
  
  setwd(dirname)
  
  #
  index_GM[[dirname]]<- runtime_GM[[dirname]] <- c()
  index_GMN[[dirname]]<- runtime_GMN[[dirname]] <- c()
  index_GMSS[[dirname]]<- runtime_GMSS[[dirname]] <- c()
  vb_it_GM[[dirname]] <-  it_GM[[dirname]] <- list()
  vb_it_GMN[[dirname]] <- it_GMN[[dirname]]  <- list()
  vb_it_GMSS[[dirname]] <- it_GMSS[[dirname]]  <- list()
  
  for(seed in 1:nrep){
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GM[[dirname]] <- c(runtime_GM[[dirname]], tmp)
    index_GM[[dirname]][[seed]] <- out$index
    it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GM[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
    # gmn
    #
    if(file.exists(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))){
      load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
      tmp <- out$total_pt
      units(tmp) <- "secs"
      runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], tmp)
      index_GMN[[dirname]][[seed]] <- out$index
      it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
      vb_it_GMN[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    }else{
      runtime_GMN[[dirname]] <- c(runtime_GMN[[dirname]], NA)
      index_GMN[[dirname]][[seed]] <- NA
      it_GMN[[dirname]][[seed]] <- NA
      vb_it_GMN[[dirname]][[seed]] <- NA
    }
    
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    tmp <- out$total_pt
    units(tmp) <- "secs"
    runtime_GMSS[[dirname]] <- c(runtime_GMSS[[dirname]], tmp)
    index_GMSS[[dirname]][[seed]] <- out$index
    it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)x$it)
    vb_it_GMSS[[dirname]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    
  }
  
}

cat("Plot runtime .. \n")
runtime_GM <- do.call("rbind", runtime_GM)
runtime_GMN <- do.call("rbind", runtime_GMN)
runtime_GMSS <- do.call("rbind", runtime_GMSS)
runtime_GM <- apply(runtime_GM, 1, function(x){c(mean(x), sd(x), length(x))})
runtime_GMN <- apply(runtime_GMN, 1, function(x){c(mean(x, na.rm = T), sd(x, na.rm = T), sum(!is.na(x)))})
runtime_GMSS <- apply(runtime_GMSS, 1, function(x){c(mean(x), sd(x), length(x))})


df <- data.table(rbind(t(runtime_GM),t(runtime_GMN),t(runtime_GMSS)), keep.rownames = T)
df$variable <- rep(c("GM","GMN","GMSS"), each = length(unique(df$rn)))
setnames(df, c(paste0("V",1:3)), c("M","SD","N"))
df[,L:=M-SD/sqrt(N)]
df[,U:=M+SD/sqrt(N)]
df[,Q := as.numeric(str_extract(rn, "(?<=q_)\\d+"))]
# save(df, file = file.path(dir_plot, "runtime_vs_q.rda"))


df[variable == "GM",variable :="GM*"]
df <- df[df$N >= ceiling(0.9 * nrep)]
ggplot(df[df$variable!="GM*" & df$Q<=100,], aes(x=Q, y=M, group=variable, color = variable,lty=variable)) +
  geom_point() +
  geom_line() +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width = 0.1) +
  scale_linetype_manual(values = 2:1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x= "number of variables (Q)", y="second",color ="",lty="")  +
  scale_color_grey(start = 0.5, end = 0.8)
ggsave(file.path(dir_plot, "runtime_vs_q_gmn_gmss.pdf"),width = 4, height = 3)


ggplot(df[df$variable!="GM*" & df$Q<=100,], aes(x=Q, y=M, group=variable, color = variable,lty=variable)) +
  geom_point() +
  geom_line() +
  geom_errorbar(mapping = aes(ymin = L, ymax=U), width = 0.1) +
  scale_linetype_manual(values = 2:1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_y_continuous(trans = "log10") +
  labs(x= "number of variables (Q)", y="second (log-scale)",color ="",lty="")  +
  scale_color_grey(start = 0.5, end = 0.8)
ggsave(file.path(dir_plot, "log_runtime_vs_q_gmn_gmss.pdf"),width = 4, height = 3)


#
cat("Load outputs .. \n")
pauc_gm <- pauc_gmn <- pauc_gmss <- list()
count <- 0
for(q in c(20, 50, 100, 200, 300, 400, 500, 600, 800, 1000)){
  count <- count + 1
  
  dirname <-  paste0("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_",Q,"_q0_3_pnorm_zeta_0.001_noise_0_beta0_0.025_randzeta")
  setwd(dirname)
  
  list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector("list",nrep)
  list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector("list",nrep)
  
  for(seed in 1:nrep){
    
    # gm
    #
    load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    bool_up <- upper.tri(net$A)
    list_GM[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GM[[seed]] <- net$A[bool_up]
    
    # gmn
    #
    if(file.exists(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))){
      load(paste0("out_GMN_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
      list_GMN[[seed]] <- out$estimates$m_delta[bool_up]
      list_true_GMN[[seed]] <- net$A[bool_up]
    }else{
      list_GMN[[seed]] <- NA
      list_true_GMN[[seed]] <- NA
    }
    
    # gmss
    #
    load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
    
  }
  
  pdf(file.path(dir_plot, paste0("q",q,"_proc_GM_GMN_GMSS.pdf")),
      width = 4, height = 4)
  
  par(pty = "s")
  plot_roc(list_GM, list_true_GM, add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1, main = "pROC curves \n (edge selection)",lty = 3)
  plot_roc(list_GMN[!is.na(list_GMN)], list_true_GMN[!is.na(list_true_GMN)], add = T, fpr_stop = 0.1, col = "#818181", spread.scale = 1, lty = 2)
  plot_roc(list_GMSS, list_true_GMSS, add = T, fpr_stop = 0.1, col = "#D4D4D4", spread.scale = 1, lty = 1)
  
  legend("bottomright",
         legend = c("GM*","GMN","GMSS"),
         lty = 3:1,
         cex = 0.7,
         col = c("black","#818181","#D4D4D4"))
  dev.off()
  
  
  pauc_gm[[count]] <- compute_pauc(list_GM, list_true_GM,  fpr_stop = 0.1, standardise = T)
  pauc_gmn[[count]] <- compute_pauc(list_GMN[!is.na(list_GMN)], list_true_GMN[!is.na(list_true_GMN)],  fpr_stop = 0.1, standardise = T)
  pauc_gmss[[count]] <- compute_pauc(list_GMSS, list_true_GMSS,  fpr_stop = 0.1, standardise = T)
  
  
  tmp <- do.call("rbind",list_GMSS_c)
  df_m_gamma <- data.table(M = apply(tmp, 2, mean),
                           SD = apply(tmp, 2, sd),
                           N = nrow(tmp),
                           ID = seq_len(ncol(tmp)))
  
  
  df_m_gamma[,L:=pmax(M-SD/sqrt(N), 0)]
  df_m_gamma[,U:=pmin(M+SD/sqrt(N),1)]
  tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = list_true_GMSS_c[[1]]!=0)
  df_m_gamma <- merge(df_m_gamma, tmp, by=c("ID"))
  
  ggplot(df_m_gamma, aes(x=ID, y=M, color=non_zero)) +
    geom_point()+
    geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
    labs(x="\n  auxiliary variables", y = "posterior inclusion probability \n",
         color = "active", title = "Auxiliary variable PPIs\n")+
    theme_bw() +
    scale_x_discrete()+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))+
    scale_color_manual(breaks = c(FALSE,TRUE),
                       labels = c("no","yes"),
                       values=c("#BABABA","black"))+
    theme(legend.position = "bottom",panel.spacing = unit(1, "lines"))+
    geom_hline(aes(yintercept = 0.5), linetype = 2) +
    theme(plot.title =  element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(dir_plot, paste0("q",q,"_PPI_GMSS_c.pdf")),
         width = 4, height = 3.5)
  
}

cat("Plot pAUCs ..\n")
df <- rbind(data.table(M = sapply(pauc_gm, mean),
                       SD = sapply(pauc_gm, sd),
                       NR = sapply(pauc_gm, length),
                       method = "GM",
                       Q = c(20, 50, 100, 200, 300, 400, 500, 600, 800, 1000)),
            data.table(M = sapply(pauc_gmn, mean),
                       SD = sapply(pauc_gmn, sd),
                       NR = sapply(pauc_gmn, length),
                       method = "GMN",
                       Q = c(20, 50, 100, 200, 300, 400, 500, 600, 800, 1000)))

df <- rbind(df, data.table(M = sapply(pauc_gmss, mean),
                           SD = sapply(pauc_gmss, sd),
                           NR = sapply(pauc_gmss, length),
                           method = "GMSS",
                           Q = c(20, 50, 100, 200, 300, 400, 500, 600, 800, 1000)))
df$L <- df$M - df$SD/sqrt(df$NR)
df$U <- df$M + df$SD/sqrt(df$NR)
#save(df, file = file.path(dir_plot, "pauc_vs_q.rda"))

df[df$method=="GM",]$ method <- "GM*"
dcast(df, Q~method, value.var = "M")
df <- df[df$NR > ceiling(0.9 * nrep), ]
ggplot(df, aes(Q, M, color = method, linetype = method)) +
  geom_point(position=position_dodge(.05)) + geom_line()+
  scale_linetype_manual(values = 3:1) +
  geom_errorbar(aes(ymin=L, ymax=U),width =.1,
                position=position_dodge(.05)) +
  theme_classic() + theme(legend.position = "bottom",
                          plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
  labs(x="number of auxiliary variables (Q)", y = "pAUC", color = "", linetype = "",title = "Edge selection performance") +
  scale_color_grey() +
  geom_vline(xintercept = 50, linetype = 3)
ggsave(file.path(dir_plot, "pauc_vs_q.pdf"),
       width = 3.5, height = 4)

ggplot(df[method!="GM*",], aes(Q, M, color = method, linetype = method)) +
  geom_point(position=position_dodge(.05)) + geom_line()+
  scale_linetype_manual(values = 2:1) +
  geom_errorbar(aes(ymin=L, ymax=U),width =.1,
                position=position_dodge(.05)) +
  theme_classic() + theme(legend.position = "bottom",
                          plot.title = element_text(size=14, face="bold", hjust = 0.5)) +
  labs(x="number of auxiliary variables (Q)", y = "pAUC", color = "", linetype = "",title = "Edge selection performance") +
  scale_color_grey(start = 0.5, end = 0.8) +
  geom_vline(xintercept = 50, linetype = 3)
ggsave(file.path(dir_plot, "pauc_vs_q_gmn_gmss.pdf"),
       width = 4, height = 4)


########################

### model selection ###

########################

cat("============== model selection criteria ============ \n")

dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
setwd(dirname)

df_model_criteria <- data.table()

for (seed in 1:nrep) {
  
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
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
),by = "seed"]

df_model_criteria <- melt(df_model_criteria, id.vars = c("v0", "seed"))
df_model_criteria_summary <- df_model_criteria[,list(M = mean(value),
                                                     SD = sd(value),
                                                     N = length(value)), by = c("v0", "variable")]
df_model_criteria_summary[,L:= M - SD/sqrt(N)]
df_model_criteria_summary[,U:= M + SD/sqrt(N)]
df_model_criteria_summary[, minM := min(M), by = "variable"]
df_model_criteria_summary[, maxM := max(M), by = "variable"]
df_model_criteria_summary[, type := "general"]
df_model_criteria_summary[!grepl("pauc",variable) & M == minM, type := "best"]
df_model_criteria_summary[variable == "pauc" & M == maxM, type := "best"]

#
ggplot(df_model_criteria_summary[!grepl("auc",variable),], aes(x = round(v0,3), y = M, color = type)) +
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U), width = 0.01)+
  theme_bw() +
  facet_wrap(factor(variable,
                    c("aic","bic","ebic"),
                    c("AIC", "BIC", "EBIC"))~.,
             scales = "free",
             nrow = 1)+
  theme(legend.position = "none") +
  labs(x = "\n spike standard deviation", y = "") +
  # scale_y_continuous(expand = c(0.01,0.01), trans="log10") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_color_manual(values = c("black","grey"))

ggsave(file.path(dir_plot, "model_selection_criteria.pdf"),
       width = 9, height = 2.5)

#
ggplot(df_model_criteria_summary[grepl("auc",variable),], aes(x = round(v0,3), y = M, color = type)) +
  geom_point() +
  geom_errorbar(aes(ymin=L, ymax=U), width = 0.01)+
  theme_bw() +
  facet_wrap(factor(variable,
                    c("aic","pauc","paucc"),
                    c("AIC", "pAUC (edge)", "pAUC (auxiliary variable)"))~.,
             scales = "free",
             nrow = 1)+
  theme(legend.position = "none") +
  # labs(x = "\n spike standard deviation", y = "") +
  labs(x = "", y = "") +
  # scale_y_continuous(trans="log10") +
  scale_color_manual(values = c("black","grey"))

ggsave(file.path(dir_plot,"model_selection_main.pdf"),
       width = 6.2, height = 2.4)

########################

### zeta sensitivity ###

########################

cat("=========== sensitivity to zeta prior ===========\n")

dirname <- "/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0"
setwd(dirname)

count <- 0
list_GMSS <- list_true_GMSS <- list_GMSS_c <- list_true_GMSS_c <- 
  pt <- index <- it <- vb_it <- vector("list", 9)

for(mn in c(25, 50, 150)){
  for(stddev in c(25, 50, 150)){
    
    count <- count + 1
    list_GMSS[[count]] <- list_true_GMSS[[count]] <- list_GMSS_c[[count]] <- list_true_GMSS_c[[count]] <-
      pt[[count]] <- index[[count]]<- it[[count]] <- vb_it[[count]]<- vector("list", nrep)
    
    for(seed in 1:nrep){
      #
      load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_",mn,"_sd_",stddev,"_gridstart_0.01_init0.rda"))
      bool_up <- upper.tri(net$A)
      #
      list_GMSS[[count]][[seed]] <- out$estimates$m_delta[bool_up]
      list_true_GMSS[[count]][[seed]] <- net$A[bool_up]
      list_GMSS_c[[count]][[seed]] <- out$estimates$m_gamma
      list_true_GMSS_c[[count]][[seed]] <- beta_true!=0
      #
      tmp <- out$total_pt
      units(tmp) <- "secs"
      pt[[count]][[seed]] <- tmp
      index[[count]][[seed]] <- out$index
      it[[count]][[seed]] <- sapply(out$full_output, function(x)x$it)
      vb_it[[count]][[seed]] <- sapply(out$full_output, function(x)max(x$vec_VB_it))
    }
  }
}

cbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#3nrep288")

pdf(file.path(dir_plot,"q50q03_proc_zeta_hyper.pdf"), width = 7 ,height = 4)

par(mfrow=c(1,2), pty="s")

# proc edge 
#
for(count in 1:length(list_GMSS)){
  
  plot_roc(list_GMSS[[count]], list_true_GMSS[[count]], col = cbPalette[(count-1)%/% 3 + 1], 
           lty = (count-1)%%3 + 1, add = ifelse(count==1, FALSE, TRUE), fpr_stop = 0.1,  spread.scale = 1, 
           main = "pROC curves \n (edge selection)")
  
}

legend("bottomright",            
       legend = c("5e-3", "1e-2", "3e-2","5e-3", "1e-2", "3e-2"),
       lty = c(rep(1,3), 1:3),
       cex = 0.5,
       col = c(cbPalette[1:3], rep("black",3)), 
       title = "sparsity (mean & sd)",ncol=2)

# proc auxiliary variable
#
for(count in 1:length(list_GMSS)){
  
  plot_roc(list_GMSS_c[[count]], list_true_GMSS_c[[count]], col = cbPalette[(count-1)%/% 3 + 1], 
           lty = (count-1)%%3 + 1, add = ifelse(count==1, FALSE, TRUE), fpr_stop = 0.1,  spread.scale = 1, 
           main = "pROC curves \n (auxiliary variable selection)")
  
}

legend("bottomright",            
       legend = c("5e-3", "1e-2", "3e-2","5e-3", "1e-2", "3e-2"),
       lty = c(rep(1,3), 1:3),
       cex = 0.5,
       col = c(cbPalette[1:3], rep("black",3)), 
       title = "sparsity (mean & sd)",ncol=2)

dev.off()

# PPI auxiliary variables
#
m_gamma <- lapply(list_GMSS_c,function(x)do.call("rbind",x))
df_m_gamma <- lapply(m_gamma, function(x)data.table(
  M=apply(x, 2, mean),
  SD=apply(x, 2, sd),
  N=nrow(x),
  ID = seq_len(ncol(x))
))
names(df_m_gamma) <- paste0(rep(c("5e-3", "1e-2", "3e-2"), each = 3), "_",
                            rep(c("5e-3", "1e-2", "3e-2"), times = 3))
df_m_gamma <- rbindlist(df_m_gamma, use.names = T,idcol = "pars")
df_m_gamma$spar <- as.numeric(gsub("_.*$","",df_m_gamma$pars))
df_m_gamma$sd <- as.numeric(gsub("^.*_","",df_m_gamma$pars))
df_m_gamma$sd <- paste0("sd = ",df_m_gamma$sd)
df_m_gamma$spar <- paste0("mean = ",df_m_gamma$spar)

df_m_gamma[,L:=pmax(M-SD/sqrt(N), 0)]
df_m_gamma[,U:=pmin(M+SD/sqrt(N),1)]
tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = beta_true!=0)
df_m_gamma <- merge(df_m_gamma, tmp, by=c("ID"))
df_m_gamma$non_zero <- factor(df_m_gamma$non_zero, c(FALSE,TRUE), c("no","yes"))

ggplot(df_m_gamma, aes(x=ID, y=M, color = non_zero)) + 
  geom_point()+
  geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
  labs(x="\n auxiliary variables", y = "posterior inclusion probability \n",
       color = "active", title = "Auxiliary variable PPIs \n")+
  theme_bw() +
  scale_x_discrete()+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  scale_color_grey(start = 0.8, end = 0.2)+
  theme(legend.position = "bottom",panel.spacing = unit(1, "lines"),
        plot.title =  element_text(hjust = 0.5, face = "bold",size = 20))+
  geom_hline(aes(yintercept = 0.5), linetype = 2) +
  facet_grid(sd~spar) 

ggsave(file.path(dir_plot,"q50q03_PPI_zeta_hyper.pdf"),
       width = 10, height = 6)

########################

###     P = 1000     ###

########################

cat("P = 1000: \n")
setwd("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_1000_q_50_q0_3_zeta_-1.55_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
pt_gm <- pt_gmss <- c()
for(seed in 1:nrep){
  load(paste0("out_GM_v2_VBECM_seed_",seed,"_zeta_mean_50_sd_150_fixed_v0_0.076_init0.rda"))
  units(out$pt) <- "hours"
  pt_gm <- c(pt_gm, out$pt)
  load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_fixed_v0_0.076_init0.rda"))
  units(out$pt) <- "hours"
  pt_gmss <- c(pt_gmss, out$pt)
}
cat(round(mean(pt_gm),2), "(",round(sd(pt_gm)/sqrt(nrep),2),")\n")
cat(round(mean(pt_gmss),2), "(",round(sd(pt_gmss)/sqrt(nrep),2),")\n")

#####################################

###   ECM & VBECM investigation   ###

#####################################
cat("ECM & VBECM investigation: \n")

setwd("/rds-d6/user/xx825/hpc-work/navigm_out_revised/simulation_n_200_p_100_q_50_q0_3_zeta_-1.52_noise_0.1_beta0_0.5_codata_TRUE_hub_TRUE_corr_0")
cat("fixed v0, s0, v1, s1 ... \n")

load("out_GMSS_zeta_mean_50_sd_150_v0_0.142_v1_100_s0_1e-06_s1_1.rda")
bool_up <- upper.tri(list_net[[1]])

cat("plot pROC curves .. \n")
pdf(file.path(dir_plot, paste0("proc_gmss_em_vbem_v0_0.142_v1_100_s0_1e-06_s1_1.pdf")), width = 8, height = 4)
par(mfrow = c(1,2), pty= "s")
# edges
#
plot_roc(lapply(list_out_vbecm, function(x)x$estimates$m_delta[bool_up]), 
         lapply(list_net, function(x)x[bool_up]), 
         add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1,
         main = "pROC curves \n (edge selection)", lty=2, ylab = "")
plot_roc(lapply(list_out_ecm, function(x)x$estimates$P1[bool_up]), 
         lapply(list_net,function(x)x[bool_up]),
         add = T, fpr_stop = 0.1,  spread.scale = 1, col = "grey")
legend("bottomright",
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2,1),
       cex = 0.7,
       col = c("black","grey"))

# auxiliary variables
#
plot_roc(lapply(list_out_vbecm, function(x)x$estimates$m_gamma), 
         lapply(list_beta_true, function(x)x!=0),
         add = F, fpr_stop = 0.1, col = "black",  spread.scale = 1,
         main = "pROC curves \n (auxiliary variable selection)", ylab = "",lty = 2)
plot_roc(lapply(list_out_ecm, function(x)x$estimates$P2), 
         lapply(list_beta_true, function(x)x!=0), 
         add = T, fpr_stop = 0.1,  spread.scale = 1, col = "grey")
legend("bottomright",
       legend = c("GMSS-VBECM","GMSS-ECM"),
       lty = c(2, 1),
       cex = 0.7,
       col = c("black","grey"))
dev.off()
#

# local mode? 
cat("local modes: \n")
cat("fixed v0, s0, v1, s1 and random initialisations for data 1 ... \n")

qf <- elbo <- c()
m_gamma_ecm <- m_gamma_vbecm <- matrix(nrow = 0, ncol = 50)
m_delta_ecm <- m_delta_vbecm <- matrix(nrow = 0, ncol = 4950)

for(seed in 1:200){
  
  load(paste0("out_GMSS_ECM_dseed_1_iseed_",seed,"_zeta_mean_50_sd_150_v0_0.142_v1_100_s0_1e-06_s1_1.rda"))
  qf <- c(qf, out$debugs$vec_ELBO_CM[length(out$debugs$vec_ELBO_CM)])
  m_delta_ecm <- rbind(m_delta_ecm, out$estimates$P1[bool_up])
  m_gamma_ecm <- rbind(m_gamma_ecm, out$estimates$P2)
  
  load(paste0("out_GMSS_VBECM_dseed_1_iseed_",seed,"_zeta_mean_50_sd_150_v0_0.142_v1_100_s0_1e-06_s1_1.rda"))
  elbo <- c(elbo, out$debugs$vec_ELBO_CM[length(out$debugs$vec_ELBO_CM)])
  m_delta_vbecm <- rbind(m_delta_vbecm, out$estimates$m_delta[bool_up])
  m_gamma_vbecm <- rbind(m_gamma_vbecm, out$estimates$m_gamma)
  
}
# save(elbo, qf, m_delta_vbecm, m_delta_ecm, m_gamma_vbecm, m_gamma_ecm, file = "~/key_quan_ref_v0_0.142_v1_100_s0_1e-06_s1_1.rda")

df <- rbind(data.table(objfun = elbo, method = "VBECM", id = 1:length(elbo)),
            data.table(objfun = qf, method = "ECM", id = 1:length(qf)))

ggplot(df,aes(objfun, fill = method)) +
  geom_histogram(alpha = 0.6) +
  theme_bw() +
  labs(x="optimal values of objective functions", fill = "", title = "Optimised objective functions\n") +
  scale_fill_grey() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 17)) 
ggsave(file.path(dir_plot, "multimod_density.pdf"), width = 4, height = 4.5)


###########################

###    misspecified     ###

###########################

vec_dirname <- c("~/simulation_rand_n_200_p_100_q_50_q0_2_zeta_-2.33_noise_0_beta0_0.4_misspecified_edge/",
                 "~/simulation_n_200_p_100_q_50_q0_3_zeta_-0.84_noise_0_beta0_-0.5_opposite/",
                 "~/simulation_n_200_p_100_q_50_q0_3_zeta_-2.33_noise_0_beta0_0.15_0_-0.5_opposite/",
                 "~/simulation_n_200_p_100_q_50_q0_3_zeta_-0.84_noise_0_beta0_-0.5_opposite/")
count = 0

for(dirname in vec_dirname){
  
  cat(dirname, "\n")
  setwd(dirname)
  
  count <- count + 1
  if(count == 1){
    label <- "misspecified_edge"
  }else if(count == 2){
    label <- "opposite_q03"
  }else if(count == 3){
    label <- "combine_q02_sparse"
  }else if(count == 4){
    label <- "opposite_q03_rightzeta"
  }
  
  list_GM <- list_GMN <- list_GMSS <- list_GMSS_c  <- vector("list",nrep)
  list_true_GM <- list_true_GMN <- list_true_GMSS <- list_true_GMSS_c <- vector("list",nrep)
  
  for(seed in 1:nrep){
    
    # gmss
    #
    if(count != 4){
      load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_50_sd_150_gridstart_0.01_init0.rda"))
    }else if(count == 4){
      load(paste0("out_GMSS_VBECM_seed_",seed,"_zeta_mean_1000_sd_150_gridstart_0.01_init0.rda"))
    }
    
    list_GMSS[[seed]] <- out$estimates$m_delta[bool_up]
    list_true_GMSS[[seed]] <- net$A[bool_up]
    list_GMSS_c[[seed]] <- out$estimates$m_gamma
    list_true_GMSS_c[[seed]] <- beta_true!=0
  }
  
  tmp <- do.call("rbind",list_GMSS_c)
  df_m_gamma <- data.table(M = apply(tmp, 2, mean),
                           SD = apply(tmp, 2, sd),
                           N = nrow(tmp),
                           ID = seq_len(ncol(tmp)))
  
  df_m_gamma[,L:=pmax(M-SD/sqrt(N), 0)]
  df_m_gamma[,U:=pmin(M+SD/sqrt(N),1)]
  tmp <- data.table(ID= 1:max(df_m_gamma$ID), non_zero = list_true_GMSS_c[[1]]!=0)
  df_m_gamma <- merge(df_m_gamma, tmp, by=c("ID"))
  
  ggplot(df_m_gamma, aes(x=ID, y=M, color=non_zero)) +
    geom_point()+
    geom_errorbar(aes(ymin=L, ymax=U), width=.2)+
    labs(x="\n  auxiliary variables", y = "posterior inclusion probability \n",
         color = "active", title = "Auxiliary variable PPIs\n")+
    theme_bw() +
    scale_x_discrete()+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))+
    scale_color_manual(breaks = c(FALSE,TRUE),
                       labels = c("no","yes"),
                       values=c("#BABABA","black"))+
    theme(legend.position = "bottom",panel.spacing = unit(1, "lines"))+
    theme(plot.title =  element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(dir_plot,  paste0(label,"_PPI_GMSS_c.pdf")), width = 4, height = 3.5)
}

