library('flowCore')
library('gplots')
library("dplyr")
library("ggplot2")
library("RColorBrewer")
#### parameters ####
setwd("./")
####load fcs file
f.name <- "Specimen_001_Sample_001_011019"  # name of the file you want to analyze, file must have extension ".FCS"

fcm <- read.FCS(paste0(f.name, '.FCS'))
fcm <- as.data.frame((exprs(fcm)))

# -------------------------------------------------------------------------
# historgam for the 610
# -------------------------------------------------------------------------
str(fcm)
rfp<-fcm$`G575-A`
bfp<-fcm$`V450-A`
# rfp<-rfp - min(rfp)
# bfp<-bfp - min(bfp)
rfp<-rfp+1
bfp<-bfp+1
bfp<-log10(bfp)
rfp<-log10(rfp)
summary(bfp)
summary(rfp)

df<-data.frame("BFP"=bfp,"RFP"=rfp)

str(df)
df2<-df

summary(df2$RFP)
summary(df2$BFP)


brewer.pal(8, "Greys")
display.brewer.pal(8, "Greys")
grey.pick<-"#BDBDBD"
green.pick<-"#1B9E77"
red.pick<-"#E7298A"
brewer.pal(8, "Purples")
purp.pick<-"#BCBDDC"
black.pick<-"#252525"

jitter <- position_jitter(width = 0, height = 0.5)


g<-ggplot(df2,(aes(x=RFP, y=BFP)))+
  geom_point(size=0.2,show.legend=F,position=jitter)+
  scale_color_manual(values = c(purp.pick,black.pick),alpha(0.1,0.9))+
  scale_y_continuous(limits = c(10,20),breaks = scales::pretty_breaks(n = 8))+
  scale_x_continuous(limits = c(11,15),breaks = scales::pretty_breaks(n = 4))
# xlim(10,14)+
# ylim(10,16)
g  + stat_density_2d(aes(fill = ..level..), geom="polygon",show.legend=F)+
  #scale_fill_gradient(low="#ABDDA4", high="#FEE08B",show.legend=F)+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="log10 (pa-mCh)", y="log10 (BFP)")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  geom_vline(xintercept=13, color=black.pick, size=1)+
  geom_hline(yintercept = 15.5, size=  1, color = black.pick)+
  # geom_text(aes(x=-6, label="Number of activated cell: 1,132", y=-6.5),  text=element_text(size=14))+
  # geom_text(aes(x=-6, label="Number of unactivated cells: 177,038", y=-6),   vjust = 1.2, text=element_text(size=14))+ 
  # theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA,size=2),
        panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))

ggsave("fcs.png", width = 10, height = 10,units = "cm" ,dpi = 360)
