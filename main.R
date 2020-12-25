rm(list = ls())
dev.off()
cat("\014")

library(MDMR)
library(stringr)
library(dplyr)
library(ADNIMERGE)
library(ggrepel)
library(caret)
library(e1071)
library(InformationValue)
library(pROC)
library(coin)




nodes_ID=read.table("D:\\yyw\\pd\\node.txt", sep = ',')
sub_table=read.csv("D:\\yyw\\hc_info.csv", sep=',') 

# Summary of table
sub_table <- droplevels(sub_table)
table(sub_table$DX)


table(sub_table$DX,sub_table$course)
sub_table %>% group_by(DX) %>% summarise(total = n(),avg_Age = mean(NMSS,na.rm=TRUE),sd_Age=sd(NMSS,na.rm=TRUE))
sub_table = sub_table[sub_table[, "DX"] != "HC",] #PD HC RBM
# Chi-squared test of sex
chisq.test(table(sub_table$course,sub_table$DX), simulate.p.value = TRUE)

# ANOVA test of age with post-hoc
anova(lm(sub_table$course~sub_table$DX))

