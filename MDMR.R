
rm(list = ls())
#dev.off()
cat("\014")

library(MDMR)
library(stringr)
library(dplyr)
library(ggrepel)
library(caret)
library(e1071)
library(InformationValue)
library(pROC)
library(coin)
#
type = 'scale_sum\\'
# invscale_mean invscale_sum scale_mean scale_sum

nodes_ID=read.table("D:\\PD_NEW\\node.txt", sep = ',')
sub_table=read.csv("D:\\PD_NEW\\all_feas.csv", sep=',') 
sub_table <- droplevels(sub_table,na.rm=TRUE)

###################################################################################################################
# use train or all
sub_proc=sub_table

#pre-define array
n=nrow(sub_proc) # number of subjects
q=nrow(nodes_ID) # number of parcels
FA_array=array(0, dim=c(q,q,n))

# define connectivity features for train
conn_name=as.vector(t(outer(nodes_ID[,2], nodes_ID[,2], paste, sep=".")))  # give names to each connectivity
idx_tria=lower.tri(matrix(NA,q,q), diag = FALSE)
conn_name_vec=conn_name[idx_tria]
Feat_all=matrix(NA,nrow =n,ncol = length(conn_name_vec))
rownames(Feat_all)=sub_proc$Subject
colnames(Feat_all)=conn_name_vec

# my data process


#
matrix_path=str_c('D:\\PD_NEW\\HC\\', type)
sub_img=list.files(matrix_path)
for (i in 1:61) {
  FA_mat=as.matrix(read.csv(str_c(matrix_path, sub_img[i]), sep = ",",header = FALSE))
  Feat_all[i,]=FA_mat[idx_tria]
  FA_array[,,i]=FA_mat
}

matrix_path=str_c('D:\\PD_NEW\\PD\\', type)
sub_img=list.files(matrix_path)
for (i in 1:55) {
  FA_mat=as.matrix(read.csv(str_c(matrix_path, sub_img[i]), sep = ",",header = FALSE))
  Feat_all[i + 61,]=FA_mat[idx_tria]
  FA_array[,,i + 61]=FA_mat
}

matrix_path=str_c('D:\\PD_NEW\\RBD\\', type)
sub_img=list.files(matrix_path)
for (i in 1:29) {
  FA_mat=as.matrix(read.csv(str_c(matrix_path, as.character(sub_img[i])), sep = ",",header = FALSE))
  print(i)
  Feat_all[i + 116,]=FA_mat[idx_tria]
  FA_array[,,i + 116]=FA_mat
}

group=sub_proc$DX
Feat_all=as.data.frame(Feat_all)
Feat_all=cbind(group,Feat_all)

# MDMR
to_do = as.data.frame(sub_proc)
#to_do = to_do[to_do[, "DX"] != "PD",]
n = nrow(to_do)
X=dplyr::select(to_do,DX,sex,age)
mdmr.p=matrix(nrow=q,ncol=ncol(X))
mdmr.F=matrix(nrow=q,ncol=1)
D=matrix(0,nrow = n,ncol = n)
mdmr.p_corrected=matrix(nrow=q,ncol=ncol(X)+1)

for (i in 1:q) {
  for (j in 1:n) {
    for (k in 1:n) {
      if (isTRUE(all.equal(FA_array[i,,j],rep(0,q)))){
        D[j,k]=0
      } else if (isTRUE(all.equal(FA_array[i,,k],rep(0,q)))){
        D[j,k]=0
      } else {
        D[j,k]=(2*(1-cor(FA_array[i,,j],FA_array[i,,k])))^0.5
      }  
    }
  }
  set.seed(24)
  mdmr.res <- summary(mdmr(X = X, D = D,perm.p = TRUE, nperm = 1000))
  print(i)
  print(mdmr.res$'Permutation p-value')
  print(mdmr.res$'Permutation p-value'[-1])
  mdmr.p[i,] <- mdmr.res$'Permutation p-value'[-1]
  mdmr.F[i,] <- mdmr.res$'Statistic'[2]
}

mdmr.p_corrected=matrix(nrow=q,ncol=ncol(X)+1)
mdmr.p_corrected[,1]=mdmr.F*(n-1)
mdmr.p_corrected[,2:(ncol(X)+1)]=sapply(1:ncol(X), function(i) p.adjust(mdmr.p[,i], method = "fdr"))

rownames(mdmr.p_corrected)=nodes_ID[,2]
colnames(mdmr.p_corrected)=c('Statistic',colnames(X))
mdmr.p_corrected=as.data.frame(mdmr.p_corrected)

# Number of significant connectomes
length(which(mdmr.p_corrected$DX<0.05))

#length(which(mdmr.old_p_corrected$HAMD<0.05))
# Index for significant connectomes
which(mdmr.p_corrected$DX<0.05)

# significant connectomes
rownames(mdmr.p_corrected)[which(mdmr.p_corrected$DX<0.05)]

##
##
##
##
##

# Effect size with parcel index=t 
# create D for parcel index=t 
t = 5
D=matrix(0,nrow = n,ncol = n)
for (j in 1:n) {
  for (k in 1:n) {
    if (isTRUE(all.equal(FA_array[t,,j],rep(0,q)))){
      D[j,k]=0
    } else if (isTRUE(all.equal(FA_array[t,,k],rep(0,q)))){
      D[j,k]=0
    } else {
      D[j,k]=(2*(1-cor(FA_array[t,,j],FA_array[t,,k])))^0.5
    }  
  }
}

G <- gower(D)
# permutation
G.list <- lapply(1:q, FUN = function(m) {
  set.seed(m)
  Y.shuf <- FA_array[t,,]
  Y.shuf[m,] <- sample(Y.shuf[m,])
  D2=matrix(0,nrow = n,ncol = n)
  for (j in 1:n) {
    for (k in 1:n) {
      if (isTRUE(all.equal(Y.shuf[,j],rep(0,q)))){
        D2[j,k]=0
      } else if (isTRUE(all.equal(Y.shuf[,k],rep(0,q)))){
        D2[j,k]=0
      } else {
        D2[j,k]=(2*(1-cor(Y.shuf[,j],Y.shuf[,k])))^0.5
      }  
    }
  }
  
  gower(D2)
})
names(G.list) <- nodes_ID[,2]
par(mar = c(5, 5, 4, 2) + 0.1)
EffectSize=delta(X = X, G = G, G.list = G.list, plot.res = T)

# mean connectivity for each seed node after MDMR
n_group=length(levels(as.factor(sub_train$DX)))
D_group=matrix(0,nrow = n_group,ncol = n_group)
FA_array_group=array(0, dim=c(q,q,n_group))
Seed_Connectome=matrix(0,nrow = n_group,ncol = q)
rownames(Seed_Connectome)=levels(as.factor(sub_train$DX))
colnames(Seed_Connectome)=nodes_ID[,2]
for (i in 1:n_group){
  FA_array_group[,,i]=apply(FA_array[,,which(sub_train$DX==levels(as.factor(sub_train$DX))[i])],1:2,mean) #calulate average FA matrix for each group 
  Seed_Connectome[i,]=FA_array_group[t,,i] # Average seed-based connectome for each group
}

# Test significant change of nodal degree for specific node
#TukeyHSD(aov(Thalamus_L ~ group, NodalDegree))
