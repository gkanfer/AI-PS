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
param$seqDepth
param$sizefactor
param$lmean
param$lOD


# -------------------------------------------------------------------------
#Parkin data can be load from:
#/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Figures/121619/Prkinscreen_plot/data_122319.RData
# Try first H1 tfeb  


sim.opts.park = RNAseq.SimOptions.2grp(ngenes = 3471, p.DE=0.2, lOD=Parkin_screen$lOD)

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                         lOD="cheung", lBaselineExpr="cheung")
sim.opts.Cheung$ngenes
sim.opts.Cheung$p.DE
length(sim.opts.Cheung$lBaselineExpr)
length(sim.opts.Cheung$lOD)
tset<-data.frame("lBaselineExpr"=sim.opts.Cheung$lBaselineExpr,"lod"=sim.opts.Cheung$lOD)
View(tset)
sim.opts.Cheung$lfc
sim.opts.Cheung$sim.seed
sim.opts.Cheung$design

# -------------------------------------------------------------------------
#'Load the data of H1 from TFEB starved expirment
#'First load the best 3 sgRNA as the baseline and take the result from edgeR 
#'for the CPM calculation and foldchange 
#'-------------------------------------------------------------------------
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Export/analysis/")
Final.working.data.rank<-readRDS("counts_sgRNAx3_h1.rds")
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
lBaselineExpr<-rowMeans(cpm(Final.working.data.rank,log=T))
lBaselineExpr<-unlist(as.numeric(lBaselineExpr))
colnames(exact.result)
lOD<-as.numeric(unlist(exact.result$table[2]))
length(lOD)
length(lBaselineExpr)
summary(lOD)
View(exact.result$table)
which(exact.result$table[1]>1.8)
ind<-which(exact.result$table[1]>1.8)
lfc<-exact.result$table[ind,1]
p.DE<-length(ind)/nrow(exact.result$table)
#how many guides defrinatioly expressed above cpmlog2 fold change of 1.8
test<-as.data.frame(cbind(row.names(exact.result),exact.result$table),stringsAsFactors = F)
p.DE<-length(which(test$logFC > 1.8))/nrow(test)
#the result is 0.0111 (1%)



sim.opts.H1 = RNAseq.SimOptions.2grp(ngenes = nrow(test), p.DE=0.01123,lOD=lOD, lBaselineExpr=lBaselineExpr,lfc = lfc)
simres = runSims(Nreps = c(3, 5, 7, 10), sim.opts=sim.opts.H1,
                 DEmethod="edgeR", nsims=20)
simres$pvalue
View(simres$pvalue)
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.4,
                      stratify.by="expr", delta=1)


summaryPower(powers)
plotPower(powers)



powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.4,
                      stratify.by="expr", target.by="lfc", delta=0.5)



