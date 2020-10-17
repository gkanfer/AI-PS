

rm(list = ls())
library(PROPER)
library("edgeR")
library("data.table")
library("dplyr")
library("svMisc")

library("data.table")
library("RColorBrewer")
library(ggplot2)
#library("plyr")
library(tidyr)
library("pheatmap")
#install.packages('VennDiagram')
library(VennDiagram)

#install.packages("ggrepel")
library(ggrepel)
library(seqLogo)
require(ggseqlogo)

#install.packages("RColorBrewer", "viridis")
library("RColorBrewer")
library( "viridis")

data("cheung")
data("bottomly")
length(param$lmean)
param$lmean
data("maqc")
setwd("./analysis")
dir()



Final.working.data.rank<-readRDS("counts_sgRNAx3_h1.rds" )
colnames(Final.working.data.rank)
m<-Final.working.data.rank
rm(Final.working.data.rank)
Final.working.data.rank<-m[,c(6:1),drop=F]
str(Final.working.data.rank)

rm(xh)
xh=new("DGEList")
xh$counts=Final.working.data.rank
xh$samples = data.frame("SampleID"=colnames(xh$counts),
                        "group"=as.factor(c(rep("Activated",3),rep("notActivated",3))),
                        "lib.size"=colSums(xh$counts),
                        "norm.factors"=rep(1,6),
                        remove.zeros=TRUE)
nrow(xh$counts)
genes.df<-as.data.frame(cbind(row.names(Final.working.data.rank),gsub("_._.*","",row.names(Final.working.data.rank),perl = T)),stringsAsFactors = F)
str(genes.df)
xh$genes=genes.df
rownames(xh$genes)<-row.names(Final.working.data.rank)
barplot(colSums(xh$counts),las=2,cex.axis = 0.8,cex.names = 0.5,ylim=c(0,10000000))
plotMDS(xh,labels = xh$samples$group,cex = 0.5)
des<-model.matrix(~xh$samples$group)

xh = estimateDisp(xh)
sqrt(xh$common.dispersion)
exact.result = exactTest(xh, pair = c("Activated", "notActivated"))
summary(exact.result)

# lBaselineExpr
#lBaselineExpr<-rowMeans(cpm(Final.working.data.rank,log=T))
#lBaselineExpr<-unlist(as.numeric(lBaselineExpr))
lBaselineExpr<-log(Final.working.data.rank)
colnames(exact.result)
elod<-estimateCommonDisp(xh)
lOD<-elod$common.dispersion
#lOD<-sqrt(xh$common.dispersion)
#lOD<-as.numeric(unlist(exact.result$table[2]))
length(lOD)
length(lBaselineExpr)
summary(lOD)
View(exact.result$table)
temp<-as.data.frame(cbind(row.names(exact.result$table),exact.result$table))

which(exact.result$table[1]>1.8)
ind<-which(exact.result$table[1]>1.8)
lfc<-exact.result$table[ind,1]
p.DE<-length(ind)/nrow(exact.result$table)
#how many guides defrinatioly expressed above cpmlog2 fold change of 1.8
test<-as.data.frame(cbind(row.names(exact.result),exact.result$table),stringsAsFactors = F)
p.DE<-length(which(test$logFC > 1.8))/nrow(test)
#the result is 0.0111 (1%)



sim.opts.H1 = RNAseq.SimOptions.2grp(ngenes = nrow(test), p.DE=p.DE,lOD=1.34, lBaselineExpr=lBaselineExpr,lfc = lfc)
simres = runSims(Nreps = c(3,5,7,9), sim.opts=sim.opts.H1,DEmethod="edgeR", nsims=100)
simres$pvalue
View(simres$pvalue)
# powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.15,
#                      stratify.by="dispersion", delta=1.69,target.by = "effectsize")
deleta<-log(2^fc.filter.decision)
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.15,
                      delta=0,target.by = "effectsize")
powers

summaryPower(powers)
plotPower(powers)


# -------------------------------------------------------------------------
# Loop for gentrating Fold change over power analysis for the 3,5,7,10

# -------------------------------------------------------------------------
ll.powerFC<-NULL
ind<-seq(0,4,0.05)
for (i in 1:length(ind)){
  powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.15,
                        delta=ind[i],target.by = "effectsize")
  x<-summaryPower(powers)
  vec<-c(x[[17]],x[[18]],x[[19]],x[[20]])
  sample.size<-c(3,5,7,9)
  ll.powerFC<-rbind(ll.powerFC,cbind(vec,rep(ind[i],4),sample.size))
  rm(x,powers,vec)
  print(i)  
}

df.powerFC<-as.data.frame(ll.powerFC,stringsAsFactors = F)
colnames(df.powerFC)<-c("Power","FC","sample_size")
str(df.powerFC)
df.powerFC$mycolor<-NA
ind<-which(df.powerFC$sample_size==3)
df.powerFC[ind,4]<-"#000000"
ind<-which(df.powerFC$sample_size==5)
df.powerFC[ind,4]<-"#A50026"
ind<-which(df.powerFC$sample_size==7)
df.powerFC[ind,4]<-"#B3E2CD"
ind<-which(df.powerFC$sample_size==9)
df.powerFC[ind,4]<-"#FDCDAC"

str(df.powerFC)

jitter <- position_jitter(width = 0, height = 0.005)

ggplot(df.powerFC,aes(x = FC,y = Power,group=sample_size))+
  #geom_histogram(binwidth = 2.5,position = 'identity')+
  geom_line(aes(color = mycolor),position=jitter,show.legend = F,size=1)+
  #scale_color_brewer(palette="Dark2")+
  coord_cartesian(xlim = c(0,4),ylim = c(0,1))+
  labs(title="", x="Log,Fold-change (effect size)", y="Power")+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  # theme(axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA,size=2),
  #       panel.background = element_blank())
  theme(axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA,size=2))
# 
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/MS/Revision/Expirment/Power_analysis/Simulation_filter3x/")
ggsave("Parkin.png", width = 10, height = 10,units = "cm" ,dpi = 360)

