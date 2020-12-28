rm(list = ls())
dev.off()
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
library(corpcor)
library(ppcor)
library(fdrtool)

df= read.table("D:\\yyw\\pd\\all_feas2_combine.csv", sep = ',',header = TRUE)

frames <- data.frame(df);#frames
frames<-subset(frames,MoCA!="NA")
D=matrix(0,40)
for (i in 1:40) {
  col_name = paste('f',i,sep='')
  pcor_stats = pcor.test(x=frames[col_name], y=frames$MoCA, z=frames[1:2],method = "spearman")
  #print(pcor_stats)
  D[i] = pcor_stats$p.value
}
p_corrected = as.vector(D)
p_corrected=as.data.frame(p_corrected)
length(which(p_corrected$p_corrected<0.05))
which(p_corrected$p_corrected<0.05)

#D = as.vector(D)
#p_corrected = fdrtool(D,statistic = 'pvalue')#$qval
#p_corrected=as.data.frame(p_corrected)
#length(which(p_corrected$p_corrected<0.05))
#which(p_corrected$p_corrected<0.05)

#p_corrected = p.adjust(D, method = "fdr")
#p_corrected=as.data.frame(p_corrected)
#length(which(p_corrected$p_corrected<0.05))
#which(p_corrected$p_corrected<0.05)

#p_corrected = p.adjust(D, method = "BH")
#p_corrected=as.data.frame(p_corrected)
#length(which(p_corrected$p_corrected<0.05))
#which(p_corrected$p_corrected<0.05)
##############

###end

#frames = as.data.frame(frames)
#pcor_stats = pcor.test(x=frames$f1, y=frames$RBDQ.HK, z=frames[1:2],method = "spearman")
#print(pcor_stats)
#D[i] = pcor_stats$p.value

