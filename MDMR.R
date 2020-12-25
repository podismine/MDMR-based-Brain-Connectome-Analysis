
rm(list = ls())
#dev.off()
cat("\014")

library(MDMR)
library(stringr)
library(dplyr)
#library(ADNIMERGE)
library(ggrepel)
library(caret)
library(e1071)
library(InformationValue)
library(pROC)
library(coin)

nodes_ID=read.table("D:\\yyw\\pd\\node.txt", sep = ',')
sub_table=read.csv("D:\\yyw\\hc_info.csv", sep=',') 
# Summary of table
sub_table <- droplevels(sub_table,na.rm=TRUE)
table(sub_table$RBDQ.HK)


table(sub_table$RBDQ.HK,sub_table$DX)
sub_table %>% group_by(RBDQ.HK) %>% summarise(total = n(),avg_Age = mean(age,na.rm=TRUE),sd_Age=sd(age,na.rm=TRUE))
# Chi-squared test of sex
chisq.test(table(sub_table$age,sub_table$RBDQ.HK), simulate.p.value = TRUE)

# ANOVA test of age with post-hoc
anova(lm(sub_table$DX~sub_table$RBDQ.HK))
TukeyHSD(x=aov(lm(sub_table$DX~sub_table$RBDQ.HK)))

###################################################################################################################
# subgroups
idx_data = createDataPartition(y=sub_table$DX, p=0.667, list=FALSE)
sub_train <- sub_table[idx_data,]
sub_test <- sub_train[-idx_data,]
sub_train = sub_table
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
rownames(Feat_all)=sub_train$Subject
colnames(Feat_all)=conn_name_vec

# my data process
matrix_path='D:\\yyw\\pd\\HC\\'
sub_img=list.files(matrix_path)
for (i in 1:61) {
  FA_mat=as.matrix(read.csv(str_c("D:\\yyw\\pd\\HC\\", sub_img[i]), sep = ",",header = FALSE))
  #print(length(FA_mat))
  
  Feat_all[i,]=FA_mat[idx_tria]
  #print(length(FA_mat))
  FA_array[,,i]=FA_mat
}

matrix_path='D:\\yyw\\pd\\PD\\'
sub_img=list.files(matrix_path)
for (i in 1:55) {
  FA_mat=as.matrix(read.csv(str_c("D:\\yyw\\pd\\PD\\", sub_img[i]), sep = ",",header = FALSE))
  #print(length(FA_mat)
  Feat_all[i + 61,]=FA_mat[idx_tria]
  FA_array[,,i + 61]=FA_mat
}

matrix_path='D:\\yyw\\pd\\RBM\\'
sub_img=list.files(matrix_path)
for (i in 1:29) {
  FA_mat=as.matrix(read.csv(str_c("D:\\yyw\\pd\\RBM\\", as.character(sub_img[i])), sep = ",",header = FALSE))
  #print(i)
  Feat_all[i + 116,]=FA_mat[idx_tria]
  FA_array[,,i + 116]=FA_mat
}

group=sub_proc$AS
Feat_all=as.data.frame(Feat_all)
Feat_all=cbind(group,Feat_all)


# create a matrix to store FA strength
NodalDegree=matrix(nrow=n ,ncol = q)
for (j in 1:n) {
  NodalDegree[j,]=rowSums(FA_array[,,j]) 
}

NodalDegree<-as.data.frame(NodalDegree)
rownames(NodalDegree)=sub_train$Subject
colnames(NodalDegree)=nodes_ID$V2


# calulate average FA strength for each group
NodalDegree_group=aggregate(NodalDegree, list(group), mean) # mean of total strength
NodalDegree_group2=aggregate(NodalDegree, list(group), sd) # Std of total strength
rownames(NodalDegree_group)=NodalDegree_group$Group.1
NodalDegree_group=NodalDegree_group[,-1]

# GLM for fiber strength
NodalDegree_rs1=matrix(nrow=q,ncol=3)
# NodalDegree_rs2=vector(mode='numeric',length = q)

for (i in 1:q) {
  # fit1.rs=summary(aov(NodalDegree[,i]~sub_train$DX+sub_train$Sex+sub_train$Age))
  # test=as.data.frame(cbind(NodalDegree[,i],sub_train$DX))
  # test %>% group_by(V2) %>% summarise(total = n(),avg = mean(V1),sd=sd(V1))
  fit1.rs=summary(lm(NodalDegree[,i]~sub_train$DX+sub_train$sex+sub_train$age))
  # NodalDegree_rs1[i,1]=fit1.rs[[1]][["t value"]][1]
  # NodalDegree_rs1[i,2]=fit1.rs[[1]][["Pr(>F)"]][1]
  NodalDegree_rs1[i,1]=fit1.rs$coefficients[2,3]
  NodalDegree_rs1[i,2]=fit1.rs$coefficients[2,4]
}
NodalDegree_rs.p_corrected=matrix(nrow=q,ncol=2)
NodalDegree_rs.p_corrected[,1]= NodalDegree_rs1[,1]
NodalDegree_rs.p_corrected[,2]= p.adjust(NodalDegree_rs1[,2], method = "fdr")
NodalDegree_rs.p_corrected=as.data.frame(NodalDegree_rs.p_corrected)
# NodalDegree_rs.p_corrected=as.data.frame(cbind(p.adjust(NodalDegree_rs1, method = "fdr"),p.adjust(NodalDegree_rs2, method = "fdr")))
rownames(NodalDegree_rs.p_corrected)=nodes_ID[,2]
colnames(NodalDegree_rs.p_corrected)=c('Statistic','p-value')


# MDMR
to_do = as.data.frame(sub_proc)
#to_do = to_do[to_do[, "DX"] != "PD",]
#to_do<-subset(to_do,RBDQ.HK!="NA")
n = nrow(to_do)
X=select(to_do,DX,sex,age)
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


library(fdrtool)
mdmr.p_corrected[,1]=mdmr.F*(n-1)
#mdmr.p_corrected[1:q,2:(ncol(X)+1)]=sapply(1:ncol(X), function(i) fdrtool(mdmr.p[,i],statistic = 'pvalue')$pval)[1:q,]
#mdmr.p_corrected[,2:(ncol(X)+1)]=sapply(1:ncol(X), function(i) p.adjust(mdmr.p[,i], method = "fdr"))
mdmr.p_corrected[1:q,2:(ncol(X)+1)]=sapply(1:ncol(X), function(i) fdrtool(mdmr.p[,i],statistic = 'pvalue')$pval)[1:q,]
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
#write.csv(mdmr.p_corrected,"F://new_result//RBDQ.HK_p.csv")

