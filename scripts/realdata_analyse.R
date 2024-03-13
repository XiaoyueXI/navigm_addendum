library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)
library(xtable)
require(snpStats)
library(LDheatmap)

setwd("~/realdata/")
fdr <- 0.2;

# load thresholds
#
load(paste0("~/Downloads/My files 2023-08-16T12_02_21/V_matrix_cedar_LYZ_region_FDR_",fdr,".RData"))

# load data
#
load(paste0("realdata_cedar_fdr",fdr,".rda"))
str(Y)
str(V)

load(paste0("realdata_cedar_fdr",fdr,"_ifr.rda"))
str(Y)
str(V)

#
print(apply(V, 2, function(x)sum(x> ppi_thres_cedar)))

#
load("~/realdata_out/out_GM_v2_VBECM_zeta_mean_50_sd_150_gridstart_0.01_init0.rda")
out_gm <- out
load("~/realdata_out/out_GM_v2_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01_ifr.rda")
out_gm_ifn <- out
load("~/realdata_out/out_GMSS_VBECM_zeta_mean_50_sd_150_gridstart_0.01_init0.rda")
out_gmss <- out
load("~/realdata_out/out_GMss_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01_ifr.rda")
out_gmss_ifn <- out

# bool_bfdr <- T
bool_bfdr <- F
if(bool_bfdr){
  ppi <- out_gm$estimates$m_delta[upper.tri(out_gm$estimates$m_delta)]
  min(ppi)
  thres_out_gm <- get_fdr_threshold(ppi, threshold_v = seq(min(ppi), 1, length.out = 1e5))
  sum(ppi > thres_out_gm)/sum(upper.tri(out_gm$estimates$m_delta))
  sum(ppi > 0.5)/sum(upper.tri(out_gm$estimates$m_delta))
  
  ppi <- out_gm_ifn$estimates$m_delta[upper.tri(out_gm_ifn$estimates$m_delta)]
  min(ppi)
  thres_out_gm_ifn <- get_fdr_threshold(ppi, threshold_v = seq(min(ppi), 1, length.out = 1e5))
  sum(ppi > thres_out_gm_ifn)/sum(upper.tri(out_gm_ifn$estimates$m_delta))
  sum(ppi > 0.5)/sum(upper.tri(out_gm_ifn$estimates$m_delta))
  
  ppi <- out_gmss$estimates$m_delta[upper.tri(out_gmss$estimates$m_delta)]
  min(ppi)
  thres_out_gmss <- get_fdr_threshold(ppi, threshold_v = seq(min(ppi), 1, length.out = 1e5))
  sum(ppi >= thres_out_gmss)/sum(upper.tri(out_gmss$estimates$m_delta))
  sum(ppi >= 0.5)/sum(upper.tri(out_gmss$estimates$m_delta))
  
  ppi <- out_gmss_ifn$estimates$m_delta[upper.tri(out_gmss_ifn$estimates$m_delta)]
  thres_out_gmss_ifn <- get_fdr_threshold(ppi, threshold_v = seq(min(ppi), 1, length.out = 1e5))
  sum(ppi >= thres_out_gmss_ifn)/sum(upper.tri(out_gmss_ifn$estimates$m_delta))
  sum(ppi >= 0.5)/sum(upper.tri(out_gmss_ifn$estimates$m_delta))
}else{
  thres_out_gm <- thres_out_gm_ifn <- thres_out_gmss <- thres_out_gmss_ifn <- 0.5
}


#
bool_up <- upper.tri(out_gm$estimates$m_delta)
cat("Identified by GM* and GMSS/Identified by GM* or GMSS: \n")
print(sum(out_gm$estimates$m_delta[bool_up] > thres_out_gm & out_gmss$estimates$m_delta[bool_up] > thres_out_gmss)/
sum(out_gm$estimates$m_delta[bool_up] > thres_out_gm | out_gmss$estimates$m_delta[bool_up] > thres_out_gmss))
cat("Identified by GM* only /Identified by GM* or GMSS: \n")
sum(out_gm$estimates$m_delta[bool_up] > thres_out_gm & out_gmss$estimates$m_delta[bool_up] < thres_out_gmss)/
  sum(out_gm$estimates$m_delta[bool_up] > thres_out_gm | out_gmss$estimates$m_delta[bool_up] > thres_out_gmss)
cat("Identified by GMSS only /Identified by GM* or GMSS: \n")
sum(out_gm$estimates$m_delta[bool_up] < thres_out_gm & out_gmss$estimates$m_delta[bool_up] > thres_out_gmss)/
  sum(out_gm$estimates$m_delta[bool_up] > thres_out_gm | out_gmss$estimates$m_delta[bool_up] > thres_out_gmss)

#
in_gm_not_in_gmss <- which(out_gm$estimates$m_delta > thres_out_gm & out_gmss$estimates$m_delta <thres_out_gmss, arr.ind = T)
in_gm_not_in_gmss <- in_gm_not_in_gmss[in_gm_not_in_gmss[,1] < in_gm_not_in_gmss[,2],]
in_gmss_not_in_gm <- which(out_gm$estimates$m_delta < thres_out_gm & out_gmss$estimates$m_delta > thres_out_gmss, arr.ind = T)
in_gmss_not_in_gm <- in_gmss_not_in_gm[in_gmss_not_in_gm[,1] < in_gmss_not_in_gm[,2],]

#
genes_cntr_rs1384 <- names(which(V[,"rs1384"]> ppi_thres_cedar))
cat("SNPs dentified by GM only controlled by rs1384? \n")
sum(rownames(V)[in_gm_not_in_gmss[,1]] %in% genes_cntr_rs1384 |
rownames(V)[in_gm_not_in_gmss[,2]] %in% genes_cntr_rs1384)/nrow(in_gm_not_in_gmss)
cat("SNPs dentified by GMSS only controlled by rs1384? \n")
sum(rownames(V)[in_gmss_not_in_gm[,1]]%in% genes_cntr_rs1384 |
rownames(V)[in_gmss_not_in_gm[,2]]%in% genes_cntr_rs1384)/nrow(in_gmss_not_in_gm)

#
sorted_changes <- sort(apply(out_gmss$estimates$m_delta, 1, function(x)sum(x>thres_out_gmss)) - apply(out_gm$estimates$m_delta, 1, function(x)sum(x>thres_out_gm)))
cat("PPIs of rs1384 regulation for genes sorted by degree changes: \n")
print(V[names(sorted_changes), out_gmss$estimates$m_gamma > thres_out_gmss])

#
sorted_changes_ifn <- sort(apply(out_gmss_ifn$estimates$m_delta, 1, function(x)sum(x>thres_out_gmss_ifn)) - apply(out_gm_ifn$estimates$m_delta, 1, function(x)sum(x>thres_out_gm_ifn)))
cat("PPIs of regulation by selected SNPs for stimulated genes sorted by degree changes: \n")
apply(V[names(sorted_changes_ifn),out_gmss_ifn$estimates$m_gamma > thres_out_gmss_ifn], 1, function(x)sum(x>thres_out_gmss_ifn))

#
cat("degrees of key genes in the literature: \n")
cat("unstimulated: \n")
print(apply(out_gm$estimates$m_delta, 1, function(x)sum(x>thres_out_gm))[c("LYZ","YEATS4","CREB1")])
print(apply(out_gmss$estimates$m_delta, 1, function(x)sum(x>thres_out_gmss))[c("LYZ","YEATS4","CREB1")])
cat("stimulated: \n")
print(apply(out_gm_ifn$estimates$m_delta, 1, function(x)sum(x>thres_out_gm_ifn))[c("LYZ","YEATS4","CREB1")])
print(apply(out_gmss_ifn$estimates$m_delta, 1, function(x)sum(x>thres_out_gmss_ifn))[c("LYZ","YEATS4","CREB1")])

# top gene in stimulated monocyte network
(V["TP53BP2",out_gmss_ifn$estimates$m_gamma > thres_out_gmss_ifn])

# venn plots for selected hotspots
#
library(ggVennDiagram)
x <- list(rs6581889 = names(which(V[,"rs6581889"] > ppi_thres_cedar)),
          rs589448 = names(which(V[,"rs589448"] > ppi_thres_cedar)),
          rs1384 = names(which(V[,"rs1384"] > ppi_thres_cedar)),
          rs10784774 = names(which(V[,"rs10784774"] > ppi_thres_cedar)))
g <- ggVennDiagram(x,
              label = "count", 
              label_color = "black",
              label_alpha = 0,
              category.names = c("      rs6581889", "rs589448", "rs1384", "rs10784774           "),
              set_size = 3,
              label_size = 3, 
              edge_lty = "dashed", 
              edge_size = 1) +
  scale_fill_gradient(low="white",high = "grey",name = "gene count \n")+
  scale_color_grey() +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5, barwidth = 6 ))
ggsave(filename = dir_plot("venn_plots.pdf"),plot = g, width = 4, height = 6)

# LD link
# European
# rs6581889,rs2168029: 0.7135
# rs6581889,rs1384: 0.783
# rs6581889,rs589448: 0.8127
# rs6581889,rs10784774: 0.7874

# rs2168029,rs1384: 0.8782
# rs2168029,rs589448: 0.5744
# rs2168029,rs10784774: 0.8891

# rs1384, rs589448: 0.6419
# rs1384, rs10784774: 0.9881

# rs589448, rs10784774: 0.645

# load outputs and make plots
#
degree_gm <- degree_gmss <- df_highlight <- df_mean <- df_median <-
  df <- active_variants <- dfes_list<- list()

count <- 0

for(label in c("", "_ifr")){
  cat("=== ",ifelse(label == "", "Unstimulated","Stimulated"), " ===\n")
  count <- count + 1
  
  cat("Load", paste0("realdata_cedar_fdr",fdr,label,".rda"), "\n")
  load(paste0("realdata_cedar_fdr",fdr,label,".rda"))
  
  cat("============ data overview ============ \n")
  print(str(V))
  print(str(Y))
  cat("Number of genes controlled by a genetic variant: \n")
  print(table(apply(V, 2, function(x)sum(x>ppi_thres_cedar))))
  print(sort(apply(V, 2, function(x)sum(x>ppi_thres_cedar))))
  ppi_thres_cedar
  
  # gm
  #
  if(label == ""){
    load(paste0("~/realdata_out/out_GM_v2_VBECM_zeta_mean_50_sd_150_gridstart_0.01",label,"_init0.rda"))
    thres <- thres_out_gm
  }else if(label == "_ifr"){
    load(paste0("~/realdata_out/out_GM_v2_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01",label,".rda"))
    thres <- thres_out_gm_ifn
  }
  
  out_gm <- out
  bool_up <- upper.tri(out_gm$estimates$m_delta)
  cat("sparsity (GM*)= ", sum(out_gm$estimates$m_delta[bool_up] > thres)/sum(bool_up), "\n")
  degree_gm[[count]] <- apply(out_gm$estimates$m_delta > thres, 2, sum)
  
  # gmss
  #
  if(label == ""){
    load(paste0("~/realdata_out/out_GMSS_VBECM_zeta_mean_50_sd_150_gridstart_0.01",label,"_init0.rda"))
    thres <- thres_out_gmss
  }else if(label == "_ifr"){
    load(paste0("~/realdata_out/out_GMSS_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01",label,".rda"))
    thres <- thres_out_gmss_ifn
  }
  out_gmss <- out
  bool_up <- upper.tri(out_gmss$estimates$m_delta)
  cat("sparsity (GMSS)= ", sum(out_gmss$estimates$m_delta[bool_up] > thres)/sum(bool_up), "\n")
  degree_gmss[[count]] <- apply(out_gmss$estimates$m_delta > thres, 2, sum)
  
  # effect sizes
  m <- out_gmss$estimates$m_gamma * out_gmss$estimates$mu_beta
  sd <- sqrt(out_gmss$estimates$m_gamma * 1/out_gmss$estimates$sig2_inv_beta + (out_gmss$estimates$m_gamma - out_gmss$estimates$m_gamma^2) * out_gmss$estimates$mu_beta^2)
  dfes <- data.table(m = m, sd = sd)
  dfes$l <- dfes$m - 1.96 * dfes$sd
  dfes$u <- dfes$m + 1.96 * dfes$sd
  dfes$name <- names(out_gmss$estimates$m_gamma)
  ggplot(dfes, aes(name, m)) +
    geom_point() +
    geom_errorbar(aes(ymin=l, ymax=u)) +
    theme_classic() + 
    theme(axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) + 
    geom_hline(yintercept = 0) +
    labs(x="",y="effect sizes")
  dfes_list[[count]] <- dfes
  
  # unique and shared
  #
  cat("GM & GMSS: \n")
  print(round(table(out_gm$estimates$m_delta[bool_up] >=thres,out_gmss$estimates$m_delta[bool_up] >=thres)))
  print(round(table(out_gm$estimates$m_delta[bool_up] >=thres,out_gmss$estimates$m_delta[bool_up] >=thres)/ (ncol(Y) * (ncol(Y) - 1)/2 ) , 6))
  
  # PPIs
  #
  cat("Selected SNPs: \n")
  print(sort(out_gmss$estimates$m_gamma))
  
  # active genetic variants 
  #
  active_variants[[count]] <- which(out_gmss$estimates$m_gamma > thres)
  
  
  # neighbours of hubs in the literature
  #
  lyzneigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["LYZ",] > thres)), "LYZ")
  yeats4neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["YEATS4",] > thres)), "YEATS4")
  creb1neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["CREB1",] > thres)), "CREB1")
  
  lyzneigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["LYZ",] > thres)), "LYZ")
  yeats4neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["YEATS4",] > thres)), "YEATS4")
  creb1neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["CREB1",] > thres)), "CREB1")
  
  # degrees
  #
  df[[count]] <- rbind(data.table(degree = degree_gm[[count]],
                                  gene = names(degree_gm[[count]]), 
                                  method = "GM"),
                       data.table(degree = degree_gmss[[count]],
                                  gene = names(degree_gmss[[count]]), 
                                  method = "GMSS"))
  
  
  df_highlight[[count]] <- df[[count]][gene%in%c("LYZ","YEATS4","CREB1"),c("gene","method","degree")]
  df_highlight[[count]]$neighbour <- F
  tmp <- data.table(gene = df_highlight[[count]]$gene,
                    method = df_highlight[[count]]$method,
                    neighbour = T)
  
  tmp <- cbind(tmp, rbind(   quantile(degree_gm[[count]][creb1neigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                             quantile(degree_gm[[count]][lyzneigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                             quantile(degree_gm[[count]][yeats4neigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                             quantile(degree_gmss[[count]][creb1neigh_gmss], probs = seq(0.25, 0.75, by = 0.25)),
                             quantile(degree_gmss[[count]][lyzneigh_gmss], probs = seq(0.25, 0.75, by = 0.25)),
                             quantile(degree_gmss[[count]][yeats4neigh_gmss], probs = seq(0.25, 0.75, by = 0.25)))
  )
  setnames(tmp, c("25%","50%","75%"), c("Ldegree","degree","Udegree"))
  df_highlight[[count]] <- rbind(df_highlight[[count]], tmp, fill = T)
  
  # network mean & median
  #
  df_mean[[count]] <- df[[count]][,mean(degree), by = "method"]
  df_median[[count]] <- df[[count]][,median(degree), by = "method"]
  
  # name fix
  df[[count]][method == "GM", method:="GM*"]
  df_highlight[[count]][method == "GM", method:="GM*"]
  df_mean[[count]][method == "GM", method:="GM*"]
  df_median[[count]][method == "GM", method:="GM*"]
}

#
dfes <- rbindlist(dfes_list, idcol = "condition")
dfes$condition <- factor(dfes$condition, c(1,2), c("unstimulated","stimulated"))
tmp <- data.table(name = colnames(V),
                  nreg = apply(V, 2, function(x)sum(x>ppi_thres_cedar)))
dfes <- merge(dfes, tmp, by = "name", all = T)
dfes$label <- paste0(dfes$name, ",", dfes$nreg)



cat("Plot regulation by rs1384 across the chromosome: \n")

# load V
load(paste0("~/realdata/realdata_cedar_fdr",fdr,".rda"))
colnames(V)

# sort(unique(annot_transcripts_repl$CHR))
cntr_rs1384 <- names(which(V[,"rs1384"] > ppi_thres_cedar))
df_snp <- annot_transcripts_repl[,c("Gene","CHR")]
df_snp <- unique(df_snp[df_snp$Gene %in% rownames(V),])

#
df_snp$cntr_by_rs1384 <- F
df_snp[df_snp$Gene %in% cntr_rs1384, ]$cntr_by_rs1384 <- T
df_snp <- as.data.table(df_snp)
df_snps <- df_snp[, list(cntr_by_rs1384 = sum(cntr_by_rs1384),
                 ncntr_by_rs1384 = sum(!cntr_by_rs1384)),
          by = "CHR"]
df_snps <- rbind(df_snps, data.table(CHR=21,
                             cntr_by_rs1384 = 0,
                             ncntr_by_rs1384 = 0))
df_snps$CHR <- as.integer(df_snps$CHR)
df_snps <- melt(df_snps, id.vars = "CHR")
ggplot(df_snps, aes(fill=factor(variable, c("cntr_by_rs1384","ncntr_by_rs1384"),c("yes","no")), y=value, x=CHR)) + 
  geom_bar(position="stack", stat="identity", width = 0.5) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(0,23), expand = c(0,0),breaks = seq(1,22,by=1)) +
  scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
  labs(x = "chromosome", y = "number of genes \n in the network", fill = "controlled by rs1384") 
ggsave(file.path(dir_plot,"genes_by_chromosome.pdf"),width = 8, height = 3)


cat("Make plots: \n")
# 
load("~/Downloads/updated_SNP_data_cedar_chr12_Xiaoyue.RData")
annot_snps_repl <- annot_chr_12[annot_chr_12$rsID %in% colnames(V),]
setnames(annot_snps_repl, c("rsID","bp"), c("snp.name","position37") )
snps_repl <- snps_cedar_chr_12

dorder <- annot_snps_repl[annot_snps_repl$snp.name %in% colnames(V),]
dorder <- dorder[order(dorder$position37),]
snps <- snps_repl[,dorder$snp.name]
gsnps <- as(snps,"SnpMatrix")
ld <- LDheatmap(gsnps, flip=TRUE, 
                name="", title=NULL, 
                LDmeasure = "r",
                add.map = T, 
                geneMapLocation = 0.01, 
                geneMapLabelX=1000)

# str(ld$LDmatrix)
dfes$name <- factor(dfes$name, colnames(gsnps@.Data))
g <- ggplot(dfes[dfes$name%in%dorder$snp.name,], aes(name, m, color = condition, shape = condition)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=l, ymax=u),position=position_dodge(0.5), width = 0.5) +
  theme_classic() + 
  scale_color_grey(start = 0.7, end = 0)+
  theme(axis.text.x = element_text (angle = 90, vjust = 1, hjust=1),
        legend.position = "bottom") + 
  geom_hline(yintercept = 0) +
  labs(x="", y="", title = "effects of genetic variants") +
  theme(plot.title =  element_text(hjust = 0.5, face = "bold"))
legend_g <- get_legend(g)
#
g <- g + theme(legend.position = "none")
ggsave(file.path(dir_plot, "effects.pdf"),g,width = 6, height = 4)
print(as_ggplot(legend_g))
ggsave(file.path(dir_plot,"effects_legend.pdf"),width = 6, height = 1)
#
pdf(file.path(dir_plot,"ldcorr.pdf"), width = 4, height = 4,  onefile = F)
grid::grid.draw(ld$LDheatmapGrob)
dev.off()

#
df <- rbindlist(df, idcol = "condition")
df$condition <- factor(df$condition, c(1,2), c("unstimulated","stimulated"))
df_mean <- rbindlist(df_mean, idcol = "condition")
df_mean$condition <- factor(df_mean$condition, c(1,2), c("unstimulated","stimulated"))
df_highlight <- rbindlist(df_highlight, idcol = "condition")
df_highlight$condition <- factor(df_highlight$condition, c(1,2), c("unstimulated","stimulated"))

#
# degree plots
#

g1 <- ggplot(df,aes(degree, color = condition, linetype = method)) +
  stat_density(geom="line",position="identity")+
  theme_bw() +
  labs(y="density", x="",color="",linetype="") +
  scale_colour_grey(start = 0.6, end = 0) + 
  scale_linetype_manual(values = 2:1) +
  theme(legend.position = "none") +
  geom_vline(data = df_mean,  aes(xintercept = V1,color = condition, linetype = method), show.legend = F) 

g2 <- ggplot(df_highlight,aes(degree, gene, 
                              color = condition,
                              linetype = method,
                              shape=factor(neighbour, c(F,T),c("itself","neighbour")))) +
  geom_point(position=position_dodge(width=.5)) + 
  geom_errorbarh(aes(xmin = Ldegree, xmax = Udegree),
                 position=position_dodge(width=.5))+
  scale_shape_manual(values = c(4,19)) +
  scale_color_grey(start = 0.7, end = 0) +
  theme_bw() +
  scale_x_continuous(limits = c(0,30), expand = c(0,0)) +
  labs(shape = "", color = "", linetype = "", y = "", x = "degrees") + 
  theme(legend.position = "bottom") +
  geom_vline(data = df_mean,  
             aes(xintercept = V1,
                 color = condition,
                 linetype = method), 
             show.legend = F)

#
vec_condition <- c("unstimulated","stimulated")
g <- list()

# the active snp
for(count in 1:2){
  #
  snps <- active_variants[[count]]
  tmp <- copy(df[df$condition == vec_condition[count],])  
  tmp$snp_controlled <- F
  for (snp in snps) {
    tmp[gene %in% names(which(V[,snp] > ppi_thres_cedar)),]$snp_controlled <- T
  }
  
  tmp <- dcast(tmp, gene + snp_controlled ~ method, value.var = "degree")
  tmp[,diff := GMSS - `GM*`]
  quantile(tmp[snp_controlled==T,]$diff)
  quantile(tmp[snp_controlled==F,]$diff)
  g[[count]] <- ggplot(tmp, aes(factor(snp_controlled, 
                                       c(FALSE, TRUE), 
                                       c("no","yes")), diff))+
    ggpubr::geom_bracket(
      xmin = "no", xmax = "yes", y.position = 8,
      label.size = 2.5,
      label = "p < 0.001"
    )+
    geom_boxplot(width = 0.2, color = ifelse(count == 1, "grey60", "black")) +
    scale_y_continuous(expand = c(0,0), limits = c(-5, 10)) + 
    theme_bw() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, lty =2) +
    labs( y ="Gene degree difference \n between the GMSS and GM* networks", 
          x = paste0("\n regulation by ", ifelse(length(snps) == 1, names(snps), paste0(names(snps)[1],"\n",paste(names(snps)[-1], collapse = ", ")))
                     ) )
  
  
  # permutation test
  #
  contr_genes <- unique(tmp[tmp$snp_controlled == T, ]$gene)
  average_contr_genes <- median((degree_gmss[[count]][contr_genes] - degree_gm[[count]][contr_genes]))
  average_perm <- c()
  set.seed(123)
  nperm <- 1e5
  for(i in 1:nperm){
    tmp <- sample(rownames(V), length(contr_genes))
    average_perm <- c(average_perm, median(degree_gmss[[count]][tmp] - degree_gm[[count]][tmp]))
  }
  
  cat("p-value is ",(sum(average_perm > average_contr_genes) + 1)/(nperm + 1),"\n")
}

legend_g2 <- get_legend(g2)
g2 <- g2 + theme(legend.position = "none")
g[[2]] <- g[[2]] + labs(y="")

pdf(paste0("realdata_cedar",label,".pdf"), width = 8, height = 10) 
g1 + g2 + g[[1]] + g[[2]] + as_ggplot(legend_g2) +
plot_layout(design = "11
                      22
                      55
                      34",
heights = c(4,3,1,3),
widths = c(4,4)) + plot_annotation(theme = theme(plot.margin = margin()))
dev.off()

cat("Make tables: \n")
## tables
#

# table of hub genes
#
hotspot_controlled1 <- names(which(V[,"rs1384"] >= ppi_thres_cedar))
hotspot_controlled2 <- c()
for(hsp in c("rs6581889","rs1384","rs589448", "rs10784774","rs10879086")){
  hotspot_controlled2 <- c(hotspot_controlled2, names(which(V[,hsp] >= ppi_thres_cedar)))
}
hotspot_controlled2 <- unique(hotspot_controlled2)


df_list <- df_gm_subset_list <- df_gmss_subset_list <- list()
for (count in 1:2) {
  df <- rbind(data.table(degree = degree_gm[[count]],
                         gene = names(degree_gm[[count]]), 
                         method = "GM"),
              data.table(degree = degree_gmss[[count]],
                         gene = names(degree_gmss[[count]]), 
                         method = "GMSS"))
  df <- dcast(df, gene ~ method, value.var = "degree")
  df[,diff := GMSS- GM]
  
  # hotspot controlled? 
  if(count == 1){
    hotspot_controlled <- hotspot_controlled1
  }else{
    hotspot_controlled <- hotspot_controlled2
  }
  
  df$hotspot<- F
  df[gene %in% hotspot_controlled, hotspot := T]
  df_list[[count]] <- df
  
  # top 10%
  df_gm_subset <- df[GM >=quantile(df$GM, prob=0.9),]
  tmp1 <- df_gm_subset$gene
  df_gm_subset[,gene2:=paste0("\\textit{",gene,"}")]
  df_gmss_subset <- df[GMSS >=quantile(df$GMSS, prob=0.9),]
  tmp2 <- df_gmss_subset$gene
  df_gmss_subset[,gene2:=paste0("\\textit{",gene,"}")]
  df_gmss_subset[hotspot==F,gene2:=paste0(gene2,"$^*$")]
  print(intersect(tmp1, tmp2))
  # in common
  print(length(intersect(df_gm_subset$gene, df_gmss_subset$gene)))
  
  # 
  df_gm_subset <- df_gm_subset[order(GM, decreasing = T),]
  df_gmss_subset <- df_gmss_subset[order(GMSS, decreasing = T),]
  
  
  df_gm_subset <- subset(df_gm_subset, select = c("gene2", "GM"))
  df_gmss_subset <- subset(df_gmss_subset, select = c("gene2", "GMSS","diff"))
  
  #
  df_gmss_subset <- rbind(df_gmss_subset, data.table(matrix(NA, 
                                                            nrow = nrow(df_gm_subset) - nrow(df_gmss_subset),
                                                            ncol = ncol(df_gmss_subset))), use.names = F)
  #
  df_gm_subset_list[[count]] <- df_gm_subset
  df_gmss_subset_list[[count]] <- df_gmss_subset
}

if(nrow(df_gm_subset_list[[1]]) > nrow(df_gm_subset_list[[2]])){
  
  tmp <- cbind(df_gm_subset_list[[2]], df_gmss_subset_list[[2]])
  tmp <- rbind(tmp,
        matrix(NA, nrow = nrow(df_gm_subset_list[[1]]) - nrow(df_gm_subset_list[[2]]),
               ncol = ncol(tmp)), use.names = F)
  tmp <- cbind(cbind(df_gm_subset_list[[1]], df_gmss_subset_list[[1]]), tmp)
  
}else if (nrow(df_gm_subset_list[[1]]) < nrow(df_gm_subset_list[[2]])){
  
  tmp <- cbind(df_gm_subset_list[[1]], df_gmss_subset_list[[1]])
  tmp <- rbind(tmp,
               matrix(NA, nrow = nrow(df_gm_subset_list[[2]]) - nrow(df_gm_subset_list[[1]]),
                      ncol = ncol(tmp)), use.names = F)
  tmp <- cbind(cbind(df_gm_subset_list[[2]], df_gmss_subset_list[[2]]), tmp)
  
}
print(xtable(tmp), include.rownames=FALSE)

# table of neighbours of hub genes in the literature
#
tmpt <- list()
count <- 0
for(label in c("","_ifr")){
  count <- count + 1
  # gm
  if(label == ""){
    load(paste0("~/realdata_out/out_GM_v2_VBECM_zeta_mean_50_sd_150_gridstart_0.01",label,"_init0.rda"))
    thres <- thres_out_gm
  }else if(label == "_ifr"){
    load(paste0("~/realdata_out/out_GM_v2_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01",label,".rda"))
    thres <- thres_out_gm_ifn
  }
  out_gm <- out
  
  # gmss
  if(label == ""){
    load(paste0("~/realdata_out/out_GMSS_VBECM_zeta_mean_50_sd_150_gridstart_0.01",label,"_init0.rda"))
    thres <- thres_out_gmss
  }else if(label == "_ifr"){
    load(paste0("~/realdata_out/out_GMSS_VBECM_seed_21_zeta_mean_50_sd_150_gridstart_0.01",label,".rda"))
    thres <- thres_out_gmss_ifn
  }
  out_gmss <- out
  
  #
  lyzneigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["LYZ",] > thres)), "LYZ")
  yeats4neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["YEATS4",] > thres)), "YEATS4")
  creb1neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta["CREB1",] > thres)), "CREB1")
  
  lyzneigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["LYZ",] > thres)), "LYZ")
  yeats4neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["YEATS4",] > thres)), "YEATS4")
  creb1neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta["CREB1",] > thres)), "CREB1")
  
  # 
  n <- max(c(length(yeats4neigh_gm),
             length(yeats4neigh_gmss),
             length(creb1neigh_gm),
             length(yeats4neigh_gmss),
             length(lyzneigh_gmss),
             length(creb1neigh_gmss)))
  
  #
  df_gm_subset <- copy(df_list[[count]])
  df_gm_subset[,gene2 := gene]
  df_gm_subset[,gene2:=paste0("\\textit{",gene2,"}")]
  
  # active hotspot controlled? 
  df_gmss_subset <- copy(df_list[[count]])
  df_gmss_subset[,gene2 := gene]
  df_gmss_subset[,gene2:=paste0("\\textit{",gene2,"}")]
  df_gmss_subset[hotspot==F,gene2:=paste0(gene2,"^*")]
  
  #
  df_gm_subset <- subset(df_gm_subset, select = c("gene2","GM","gene"))
  df_gmss_subset <- subset(df_gmss_subset, select = c("gene2","GMSS","diff","gene"))
  
  tmpt[[count]] <- cbind(
    rbind(df_gm_subset[gene%in% lyzneigh_gm, c("gene2", "GM")][order(-GM),], data.table(matrix(NA, nrow = n-length(lyzneigh_gm), ncol=2)),use.names =F),
    rbind(df_gmss_subset[gene%in% lyzneigh_gmss, c("gene2", "GMSS", "diff")][order(-GMSS),], data.table(matrix(NA, nrow = n-length(lyzneigh_gmss), ncol=3)),use.names =F),
    
    rbind(df_gm_subset[gene%in% yeats4neigh_gm, c("gene2", "GM")][order(-GM),], data.table(matrix(NA, nrow = n-length(yeats4neigh_gm), ncol=2)),use.names =F),
    rbind(df_gmss_subset[gene%in% yeats4neigh_gmss, c("gene2", "GMSS", "diff")][order(-GMSS),], data.table(matrix(NA, nrow = n-length(yeats4neigh_gmss), ncol=3)),use.names =F),
    
    rbind(df_gm_subset[gene%in% creb1neigh_gm, c("gene2", "GM")][order(-GM),], data.table(matrix(NA, nrow = n-length(creb1neigh_gm), ncol=2)),use.names =F),
    rbind(df_gmss_subset[gene%in% creb1neigh_gmss, c("gene2", "GMSS", "diff")][order(-GMSS),], data.table(matrix(NA, nrow = n-length(creb1neigh_gmss), ncol=3)),use.names =F)
  )
}


print(xtable(tmpt[[1]]), include.rownames=FALSE)
print(xtable(tmpt[[2]]), include.rownames=FALSE)
