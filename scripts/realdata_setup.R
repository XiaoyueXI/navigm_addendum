library(data.table)

# load filtered auxiliary variable matrix
#
fdr <- 0.2
load(paste0("~/Downloads/My files 2023-08-16T12_02_21/V_matrix_cedar_LYZ_region_FDR_",fdr,".RData"))



# duplicates
# take average, conservative (if one is insignificant and another is significant, the average leads to insignificant)
#
dup_genes <- names(which(table(rownames(V_cedar)) > 1))
tmp <- V_cedar[rownames(V_cedar)%in% dup_genes,]
# apply(tmp, 1, function(x)which(x>=ppi_thres_cedar))
#
tmp_m <- aggregate(tmp, by=list(rownames(tmp)), FUN=mean)
rownames(tmp_m) <- tmp_m[,1]
tmp_m <- tmp_m[,-1]
stopifnot(all(colnames(tmp_m) == colnames(V_cedar)))
V <- as.matrix(rbind(V_cedar[!rownames(V_cedar)%in% dup_genes,], tmp_m))
print(str(V))


# filter after coping with duplicates
#
id2 <- which( apply(V, 2, function(x)sum(x>ppi_thres_cedar)) != 0 )
id1 <- which( apply(V, 1, function(x)sum(x>ppi_thres_cedar)) != 0 )
V <- V[id1, id2]
print(str(V))


# load expression data
#
load("~/Downloads/expression_data_monocytes_and_bcells_with_eQTL_PPI_LYZ_region.RData")
Y <- expr_unstim[, rownames(V)]
print(str(Y))


# check
#
stopifnot(all(rownames(V) == colnames(Y)))

#
save(Y, V, file = paste0('~/realdata/realdata_cedar_fdr',fdr,'.rda'))
