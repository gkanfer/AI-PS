library(EBImage)
library(e1071)
library(data.table)
setwd("C:/outproc")
if (file.exists("Binary1.tif"))  
  file.remove("Binary1.tif")
#loading model
setwd("./model")
model<-readRDS( "model_Parkin.rds")
################################################################################################################
# Loading image
################################################################################################################
frame<-readImage("C:/outproc/color.tif",type = "tif")
GFP<-frame[1:512,1:512,1]
dapi<-frame[1:512,1:512,2]
# minmax normalization
minDapi<-min(as.vector(dapi))
maxDapi<-max(as.vector(dapi))
minGFP<-min(as.vector(GFP))
maxGFP<-max(as.vector(GFP))
dapin<-normalize(dapi, ft=c(0,1),c(minDapi,maxDapi))
GFPn<-normalize(GFP, ft=c(0,1) ,c(minGFP,maxGFP))
################################################################################################################
#nuclar segmentation
################################################################################################################
dapi_normal<- dapin*2
nmask2 = thresh(dapi_normal, 16, 16,0.04)
##filling holes and remove the small obgect by open
mk3 = makeBrush(3, shape= "diamond")
nmask2 = opening(nmask2, mk3) 
nmask2 = fillHull(nmask2) 
#extraction
nmask2 = bwlabel(nmask2)  #binery to object
nf = computeFeatures.shape(nmask2)
#removing small obgects: remove all the object smaller than 50 pixel beased on the perimeter which is column 2
nr = which(nf[,2] < 30) 
nmask= rmObjects(nmask2, nr) 
nmask = watershed( distmap(nmask), 2 ) ##method to spreat two adgcent objects
#computes the geometric features of the nuclei and #prints them-gives you some data about the nucluas
nf = computeFeatures.shape(nmask)
#removing small obgects: remove all the object smaller than 50 pixel beased on the perimeter which is column 2
nr = which(nf[,2] < 50) 
nseg = rmObjects(nmask, nr) 
nn = max(nseg) 
colorMode(nseg)<-Grayscale
if ((length(nf[,1])>2000)|(length(nf[,1])<1))
{
  xr<-2:500
  xr<-sample(xr,3) 
  yr<-2:500
  yr<-sample(yr,3)
  global<-opening( dapi_normal<0.0014)
  dev.off()
  tiff("blank.tif",512,512)
  dev.off()
  blank<-readImage("blank.tif")
  blank<- Image(blank,c(512,512),"Grayscale")
  seg_point = paintObjects(blank, global, opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  # #dsiplay(seg_point/0.1, method="raster")
  
  seg_point<-thresh(seg_point) ##makes it binery
  #dsiplay(seg_point,method="raster")
  
  setwd("C:/outproc/")
  writeImage(x = seg_point, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
  # write(x,"C:/outproc/firstif")
  # 
  
}
######################################################################################
########Cell border detection
######################################################################################
cell_normal<- GFPn*2
smooth<-makeBrush(size=9,shape = "disc") 
cell_normal<-filter2(cell_normal,smooth, boundary = c("circular", "replicate"))
thr_cell<-thresh(cell_normal, w=10, h=10, offset=0.0001)
colorMode(thr_cell)<-Grayscale
cmask_ther_cell = paintObjects(thr_cell,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
cmask = opening(thr_cell, kern=makeBrush(7,shape="disc"))
cmask_ther_cell = paintObjects(cmask,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)

open2<-opening(cell_normal>9)
colorMode(open2)<-Grayscale
open_2 = paintObjects(open2,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
combine<-cmask
combine[open2 > cmask]<-open2[open2 > cmask]
combine[nmask2 > combine]<-nmask2[nmask2 > combine]
open_2 = paintObjects(combine,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
csegpink = propagate(cell_normal, nseg, lambda=1.0e-2, mask=cmask)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
csegpink = propagate(cell_normal, nseg, lambda=1.0e-2, mask=combine)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
##computes the geometric features of the cells outline
cf = computeFeatures.shape(csegpink)
if (length(cf[,1])<1)
{
  xr<-2:500
  xr<-sample(xr,3) 
  yr<-2:500
  yr<-sample(yr,3)
  global<-opening( dapi_normal<0.0014)
  dev.off()
  tiff("blank.tif",512,512)
  #dsiplay(global,method="raster")
  #celltext = text(x= xr, y= yr , labels="x", col="yellow", cex = 0.4)
  dev.off()
  
  blank<-readImage("blank.tif")
  
  blank<- Image(blank,c(512,512),"Grayscale")
  write(x,"C:/outproc/test5")  
  seg_point = paintObjects(blank, global, opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  # display(seg_point/0.1, method="raster")
  
  seg_point<-thresh(seg_point) ##makes it binery
  #dsiplay(seg_point,method="raster")
  
  setwd("C:/outproc/")
  writeImage(x = seg_point, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
  # # write(x,"C:/outproc/secondtif")
}   
##removing small obgects: remove all the object smaller (perimeter) than 50 pixel beased on the perimeter which is column 2
cr = which(cf[,2] < 200) 
csegpink = rmObjects(csegpink, cr) 
##removing small obgects: remove all the object smaller (perimeter) than 250 pixel beased on the perimeter which is column 2
cf = computeFeatures.shape(csegpink)
cfErea<-data.frame(cf[,1])
cfErea$num<-row.names(cfErea)
ci = which(cf[,1] > 12000) 
csegpink = rmObjects(csegpink, ci) 
csegpink   = propagate(csegpink, nseg, lambda=1.0e-2, mask=combine)
csegpink <- fillHull(csegpink)
colorMode(csegpink)<-Grayscale
#remove tuching objects at image edges
dims<-dim(csegpink)
border<-c(csegpink[1:dims[1],1], csegpink[1:dims[1],dims[2]], csegpink[1,1:dims[2]], csegpink[dims[1],1:dims[2]])
ids<-unique(border[which(border !=0)])
csegpink<-rmObjects(csegpink, ids)
nseg<-rmObjects(nseg, ids)
xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]

#Fetures extraction
table_test_pink_shape = computeFeatures.shape(csegpink,GFP)
table_test_pink_moment = computeFeatures.moment(csegpink,GFP)
table_test_pink_basic = computeFeatures.basic(csegpink,GFP)
table_test_pink<-data.frame(cbind(table_test_pink_basic,table_test_pink_moment,table_test_pink_shape))
rownameTable<-row.names(table_test_pink)
table_test_pink<-data.frame(cbind(rownameTable,table_test_pink))
Ts.mix<-table_test_pink[,2:12]
rowNameTable<-table_test_pink[,1]
######################################################################################
########Phnotype prediction
######################################################################################
y.pred<-predict(model,Ts.mix, decision.values = T)
#length(y.pred)
d=attr(y.pred,"decision.values")
new.y.pred=rep("p",length(y.pred))
NewCutoff=0
new.y.pred[d<NewCutoff]="N"
#length(new.y.pred)
d<-round(d,1)
Ts.mix$pred<-as.array(new.y.pred)
Ts.mix<-Ts.mix[1:length(table_test_pink[,1]),]
Ts.mix$rowNameTable<-rowNameTable
ir = which(Ts.mix$pred %in% "p") 
csegpink = rmObjects(csegpink, ir) 
pr = which(Ts.mix$pred %in% "N")
mseg<- rmObjects(nseg, nr)  
nr = which(Ts.mix$pred %in% "p") 
nseg = rmObjects(nseg, nr) 
if (length(nr) == length(Ts.mix$pred)){
  xr<-2:500
  xr<-sample(xr,3) 
  yr<-2:500
  yr<-sample(yr,3)
  global<-opening( dapi_normal<0.0014)
  dev.off()
  tiff("blank.tif",512,512)
  # display(global,method="raster")
  #celltext = text(x= xr, y= yr , labels="x", col="yellow", cex = 0.4)
  dev.off()
  
  blank<-readImage("blank.tif")
  
  blank<- Image(blank,c(512,512),"Grayscale")
  write(x,"C:/outproc/test5")  
  seg_point = paintObjects(blank, global, opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
  # display(seg_point/0.1, method="raster")
  
  seg_point<-thresh(seg_point) ##makes it binery
  #dsiplay(seg_point,method="raster")
  
  setwd("C:/outproc/")
  writeImage(x = seg_point, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
  
  
  
}
mk9 = makeBrush(9, shape= "diamond")
nnseg<-erode(nseg,mk9)
nseg_bin<-thresh(nnseg) ##makes it binery
setwd("C:/outproc/")
writeImage(x = nseg_bin, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
closeAllConnections()
