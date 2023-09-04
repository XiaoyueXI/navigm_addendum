library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)
library(xtable)

#
setwd('~/realdata/')

# load processed outputs using CEDAR and model inputs
#
fdr <- 0.2
load(paste0("~/Downloads/My files 2023-08-16T12_02_21/V_matrix_cedar_LYZ_region_FDR_",fdr,".RData"))
load(paste0('realdata_cedar_fdr',fdr,'.rda'))

# data overview
#
str(V)
str(Y)
table(apply(V, 2, function(x)sum(x>ppi_thres_cedar)))
sort(apply(V, 2, function(x)sum(x>ppi_thres_cedar)))
V['LYZ',]>ppi_thres_cedar

# gm
#
load(paste0('out_GM_v2_VBECM_zeta_mean_50_sd_150.rda'))
out_gm <- out
bool_up <- upper.tri(out_gm$estimates$m_delta)
cat('sparsity (GM*)= ', sum(out_gm$estimates$m_delta[bool_up] > 0.5)/sum(bool_up), '\n')
degree_gm <- apply(out_gm$estimates$m_delta > 0.5, 2, sum)

# gmss
#
load(paste0('out_GMSS_VBECM_zeta_mean_50_sd_150.rda'))
out_gmss <- out
bool_up <- upper.tri(out_gmss$estimates$m_delta)
cat('sparsity (GMSS)= ', sum(out_gmss$estimates$m_delta[bool_up] > 0.5)/sum(bool_up), '\n')
degree_gmss <- apply(out_gmss$estimates$m_delta > 0.5, 2, sum)

# unique and shared
#
round(table(out_gm$estimates$m_delta[bool_up] >=0.5,out_gmss$estimates$m_delta[bool_up] >=0.5))
round(table(out_gm$estimates$m_delta[bool_up] >=0.5,out_gmss$estimates$m_delta[bool_up] >=0.5)/ (ncol(Y) * (ncol(Y) - 1)/2 ) , 6)

# PPIs
#
sort(out_gmss$estimates$m_gamma)

# active genetic variants 
#
active_variants <- which(out_gmss$estimates$m_gamma >0.5)
out_gmss$estimates$m_gamma[active_variants]

# degrees
#
df <- rbind(data.table(degree = degree_gm,
                       gene = names(degree_gm), 
                       method = 'GM'),
            data.table(degree = degree_gmss,
                       gene = names(degree_gmss), 
                       method = 'GMSS'))

# neighbours of hubs in the literature
#
lyzneigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['LYZ',] > 0.5)), 'LYZ')
yeats4neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['YEATS4',] > 0.5)), 'YEATS4')
creb1neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['CREB1',] > 0.5)), 'CREB1')

lyzneigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['LYZ',] > 0.5)), 'LYZ')
yeats4neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['YEATS4',] > 0.5)), 'YEATS4')
creb1neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['CREB1',] > 0.5)), 'CREB1')

df_highlight <- df[gene%in%c('LYZ','YEATS4','CREB1'),c('gene','method','degree')]
df_highlight$neighbour <- F
tmp <- data.table(gene = df_highlight$gene,
                  method = df_highlight$method,
                  neighbour = T)

tmp <- cbind(tmp, rbind(   quantile(degree_gm[creb1neigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                           quantile(degree_gm[lyzneigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                           quantile(degree_gm[yeats4neigh_gm], probs = seq(0.25, 0.75, by = 0.25)),
                           quantile(degree_gmss[creb1neigh_gmss], probs = seq(0.25, 0.75, by = 0.25)),
                           quantile(degree_gmss[lyzneigh_gmss], probs = seq(0.25, 0.75, by = 0.25)),
                           quantile(degree_gmss[yeats4neigh_gmss], probs = seq(0.25, 0.75, by = 0.25)))
)
setnames(tmp, c("25%","50%","75%"), c('Ldegree','degree','Udegree'))
df_highlight <- rbind(df_highlight, tmp, fill = T)

# network mean & median
#
df_mean <- df[,mean(degree), by = 'method']
# df[,median(degree), by = 'method']

# name fix
df[method == 'GM', method:='GM*']
df_highlight[method == 'GM', method:='GM*']
df_mean[method == 'GM', method:='GM*']

#
g1 <- ggplot(df,aes(degree, color = method, linetype = method)) +
  stat_density(geom="line",position="identity")+
  theme_bw() +
  labs(y='density', x='',color='',linetype='') +
  scale_colour_grey() + 
  scale_linetype_manual(values = 2:1) +
  theme(legend.position = 'bottom') +
  geom_vline(data = df_mean,  aes(xintercept = V1,color = method, linetype = method), show.legend = F) +
  geom_vline(xintercept = 4,color = 'black', linetype = 3, show.legend = F)

g2 <- ggplot(df_highlight,aes(degree, gene, 
                              color = method, 
                              shape=factor(neighbour, c(F,T),c('itself','neighbour')))) +
  geom_point(position=position_dodge(width=.5)) + 
  geom_errorbarh(aes(xmin = Ldegree, xmax = Udegree),
                 position=position_dodge(width=.5))+
  scale_shape_manual(values = c(4,19)) +
  scale_color_grey() +
  scale_linetype_manual(values = c(2:1)) +
  theme_bw() +
  scale_x_continuous(limits = c(0,30), expand = c(0,0)) +
  labs(shape = '',color='',y='',x='degrees') + 
  theme(legend.position = 'bottom') +
  geom_vline(data = df_mean,  aes(xintercept = V1,color = method, linetype = method), show.legend = F)+
  geom_vline(xintercept = 4,color = 'black', linetype = 3, show.legend = F)

# pdf("degree_distribution.pdf", width = 6, height = 5) 
# g1 + g2 + plot_layout(ncol=1,heights = c(2.25,1.25), guides='collect' )&
#   theme(legend.position='bottom')
# dev.off()
# 
# g1 + theme(legend.position='bottom') + labs(x='degrees')
# ggsave("degree_distribution_top.pdf", width = 6, height = 2.5) 
# 
# g2 + theme(legend.position='bottom')
# ggsave("degree_distribution_bottom.pdf", width = 6, height = 1.5) 


# genes regulated by each of hotspot genetic variants
# and number of active genetic variants
#

# the active snp
#
snp <- "rs1384"
tmp <- copy(df)  
tmp$snp_controlled <- F
tmp[gene %in% names(which(V[,snp] > ppi_thres_cedar)),
    snp_controlled := T]
tmp <- dcast(tmp, gene + snp_controlled ~ method, value.var = 'degree')
tmp[,diff := GMSS - `GM*`]
quantile(tmp[snp_controlled==T,]$diff)
quantile(tmp[snp_controlled==F,]$diff)
g <- ggplot(tmp, aes(factor(snp_controlled, 
                            c(FALSE, TRUE), 
                            c('no','yes')), diff))+
  ggpubr::geom_bracket(
    xmin = "no", xmax = "yes", y.position = 8,
    label.size = 2.5,
    label = "p < 0.001"
  )+
  geom_boxplot(width = 0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(-5, 10)) + 
  theme_bw() +
  theme(legend.position = 'bottom') +
  geom_hline(yintercept = 0, lty =2) +
  labs( y ='Gene degree difference \n between the GMSS and GM* networks', x = paste0('\n regulation by ', snp) )


# permutation test
#
contr_genes <- unique(tmp[tmp$snp_controlled == T, ]$gene)
average_contr_genes <- median((degree_gmss[contr_genes] - degree_gm[contr_genes]))
average_perm <- c()
set.seed(123)
nperm <- 1e5
for(i in 1:nperm){
  tmp <- sample(rownames(V), length(contr_genes))
  average_perm <- c(average_perm, median(degree_gmss[tmp] - degree_gm[tmp]))
}

cat('p-value is ',(sum(average_perm > average_contr_genes) + 1)/(nperm + 1),'\n')

contr_genes_id <- which(rownames(V) %in% contr_genes)
print(exactRankTests::perm.test(degree_gmss[contr_genes_id] - degree_gm[contr_genes_id],
                                degree_gmss - degree_gm))

# uncontr_genes_id <- which(!genes %in% contr_genes)
# print(exactRankTests::perm.test(degree_gmss[contr_genes_id] - degree_gm[contr_genes_id],
#                           degree_gmss[uncontr_genes_id] - degree_gm[uncontr_genes_id]))
# ggsave(g, filename = "realdata_GMSS_GM_degree_diff_by_regulation.pdf",width = 4, height = 3.5)


pdf("realdata_cedar.pdf", width = 8, height = 5) 
g1 + g2 + g + plot_layout(design = "1#
                                    13
                                    23",
                          heights = c(1.25,1,1.25),
                          widths = c(2.25,1.25) )
dev.off()


## tables
#

# table of hub genes
#
df <- rbind(data.table(degree = degree_gm,
                       gene = names(degree_gm), 
                       method = 'GM'),
            data.table(degree = degree_gmss,
                       gene = names(degree_gmss), 
                       method = 'GMSS'))
df <- dcast(df, gene ~ method, value.var = 'degree')
df[,diff := GMSS- GM]

# hotspot controlled? 
hotspot <-  'rs1384'
hotspot_controlled <- names(which(V[,hotspot] >= ppi_thres_cedar))
df$hotspot<- F
df[gene %in% hotspot_controlled, hotspot := T]

# top 10%
df_gm_subset <- df[GM >=quantile(df$GM, prob=0.9),]
df_gm_subset[,gene2:=paste0('\\textit{',gene,'}')]
df_gmss_subset <- df[GMSS >=quantile(df$GMSS, prob=0.9),]
df_gmss_subset[,gene2:=paste0('\\textit{',gene,'}')]
df_gmss_subset[hotspot==F,gene2:=paste0(gene2,'$^*$')]

# in common
length(intersect(df_gm_subset$gene, df_gmss_subset$gene))

# 
df_gm_subset <- df_gm_subset[order(GM, decreasing = T),]
df_gmss_subset <- df_gmss_subset[order(GMSS, decreasing = T),]

# nrow(df_gm_subset)
# nrow(df_gmss_subset)

df_gm_subset <- subset(df_gm_subset, select = c('gene2', 'GM'))
df_gmss_subset <- subset(df_gmss_subset, select = c('gene2', 'GMSS','diff'))

#
df_gmss_subset <- rbind(df_gmss_subset, data.table(matrix(NA, 
                                                          nrow = nrow(df_gm_subset) - nrow(df_gmss_subset),
                                                          ncol = ncol(df_gmss_subset))), use.names = F)
#
print(xtable(cbind(df_gm_subset, df_gmss_subset)), include.rownames=FALSE)

# table of neighbours of hub genes in the literature
#
lyzneigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['LYZ',] > 0.5)), 'LYZ')
yeats4neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['YEATS4',] > 0.5)), 'YEATS4')
creb1neigh_gmss <- setdiff(names(which(out_gmss$estimates$m_delta['CREB1',] > 0.5)), 'CREB1')

lyzneigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['LYZ',] > 0.5)), 'LYZ')
yeats4neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['YEATS4',] > 0.5)), 'YEATS4')
creb1neigh_gm <- setdiff(names(which(out_gm$estimates$m_delta['CREB1',] > 0.5)), 'CREB1')

# 
n <- max(c(length(yeats4neigh_gm),
           length(yeats4neigh_gmss),
           length(creb1neigh_gm),
           length(yeats4neigh_gmss),
           length(lyzneigh_gmss),
           length(creb1neigh_gmss)))

#
df_gm_subset <- copy(df)
df_gm_subset[,gene2 := gene]
df_gm_subset[,gene2:=paste0('\\textit{',gene2,'}')]

# active hotspot controlled? 
df_gmss_subset <- copy(df)
df_gmss_subset[,gene2 := gene]
df_gmss_subset[,gene2:=paste0('\\textit{',gene2,'}')]
df_gmss_subset[hotspot==F,gene2:=paste0(gene2,'^*')]

#
df_gm_subset <- subset(df_gm_subset, select = c('gene2','GM','gene'))
df_gmss_subset <- subset(df_gmss_subset, select = c('gene2','GMSS','diff','gene'))

tmpt <- cbind(
  rbind(df_gm_subset[gene%in% lyzneigh_gm, c('gene2', 'GM')][order(-GM),], data.table(matrix(NA, nrow = n-length(lyzneigh_gm), ncol=2)),use.names =F),
  rbind(df_gmss_subset[gene%in% lyzneigh_gmss, c('gene2', 'GMSS', 'diff')][order(-GMSS),], data.table(matrix(NA, nrow = n-length(lyzneigh_gmss), ncol=3)),use.names =F),
  
  rbind(df_gm_subset[gene%in% yeats4neigh_gm, c('gene2', 'GM')][order(-GM),], data.table(matrix(NA, nrow = n-length(yeats4neigh_gm), ncol=2)),use.names =F),
  rbind(df_gmss_subset[gene%in% yeats4neigh_gmss, c('gene2', 'GMSS', 'diff')][order(-GMSS),], data.table(matrix(NA, nrow = n-length(yeats4neigh_gmss), ncol=3)),use.names =F),
  
  rbind(df_gm_subset[gene%in% creb1neigh_gm, c('gene2', 'GM')][order(-GM),], data.table(matrix(NA, nrow = n-length(creb1neigh_gm), ncol=2)),use.names =F),
  rbind(df_gmss_subset[gene%in% creb1neigh_gmss, c('gene2', 'GMSS', 'diff')][order(-GMSS),], data.table(matrix(NA, nrow = n-length(creb1neigh_gmss), ncol=3)),use.names =F)
)

#
print(xtable(tmpt), include.rownames=FALSE)

