################################################################################################################
# While loop "on the fly" AI-PS 
################################################################################################################
library(tictoc)
library(EBImage)
library(data.table)
library(outliers)
library(parallel)
library("gtools")
library(reticulate)
use_condaenv("tf-keras")
library(keras)
library(tensorflow)
library("magick")
k_clear_session()
setwd("./model")
model<-load_model_hdf5( "model_TFEB.h5" ,compile = T)
################################################################################################################
#function excuate error
################################################################################################################
next_fov<-function(){
  while (!sync.txt == 2) { 
    try(
      {
        
        #closeAllConnections()
        #setwd("c:/outproc/")  
        fileConn<-file("c:/outproc/sync.txt")
        writeLines(c("2"), fileConn)
        close(fileConn)
        print("sync.txt == 2")
        break
      }  
      , silent = TRUE
    )
  }
  
}
while (1) {
  try(
    {
      StartTime<-Sys.time()
      Sys.sleep(0.25)
      closeAllConnections()
      setwd("c:/outproc/")
      if (file.exists("Binary1.tif"))  
        file.remove("Binary1.tif") 
      sync.txt = 0
      while (1)
      {
        fileConn<-try(file("c:/outproc/sync.txt"))
        if (length(fileConn) != 0)
        {
          #sync.txt <- fread("sync.txt")
          sync.txt <- try(readLines("c:/outproc/sync.txt", 1))
          close(fileConn)
        }
        if (sync.txt == 1)
        {
          #print("sync.txt == 1")
          break
        }
      }
      # load captured image 
      frame<-readImage("C:/outproc/color.tif",type = "tif",all = T)
      dim.y<-dim(frame)[1]
      dim.X<-dim(frame)[2]
      GFP<-frame[1:dim.y,1:dim.X,1]
      dapi<-frame[1:dim.y,1:dim.X,2]
      colorMode(GFP)<-"Grayscale"
      ################################################################################################################
      #Nuclar segmentation
      ################################################################################################################
      dapi_normal<- dapi*5
      nmask2 = thresh(dapi_normal, 40, 40,0.009)
      mk3 = makeBrush(7, shape= "diamond")
      nmask2 = opening(nmask2, mk3) 
      nmask2 = fillHull(nmask2) 
      nseg = bwlabel(nmask2)  #binery to object
      seg_pink = paintObjects(nseg,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
      display(seg_pink,"raster")
      rm(nmask2)
      chackpoint<-computeFeatures.shape(nseg)
      nmask<-nseg
      if (length(chackpoint) < 3 ){
        next_fov()
      }
      if (is.null(nrow(chackpoint))){
        next_fov()
      }
      if (nrow(chackpoint) > 2000 ){
        next_fov()
      }
      nf = computeFeatures.shape(nmask)
      nr = which(nf[,2] < 120) 
      nseg = rmObjects(nmask, nr) 
      nn = max(nseg) 
      colorMode(nseg)<-Grayscale
      chackpoint<-computeFeatures.shape(nseg)
      if (length(chackpoint) < 3 ){
        next_fov()
      }
      if (is.null(nrow(chackpoint))){
        next_fov()
      }
      if (nrow(chackpoint) > 2000 ){
        next_fov()
      }
      #remove outliers 
      int.dapi<-computeFeatures.basic(nseg,dapi)
      y<-which(scores(int.dapi[,1], type="z",prob = 0.95))
      tp<-as.numeric(attr(y,"names"))
      if (length(tp) < 1 ){
        next_fov()
      }
      nseg<-rmObjects(nseg,tp)
      chackpoint<-computeFeatures.shape(nseg)
      df<-as.data.frame(chackpoint)
      rm(chackpoint,nmask)
      xy<-computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
      if (length(xy) < 3 ){
        next_fov()
      }
      if (is.null(nrow(xy))){
        next_fov()
      }
      if (nrow(xy) > 500 ){
        next_fov()
      }
      gsegg=nseg
      ######################################################################################
      #Cell border detection
      ######################################################################################
      cell_normal<- GFP*40
      smooth<-makeBrush(19,shape = "disc") 
      cell_normal<-filter2(cell_normal,smooth, boundary = c("circular", "replicate"))
      thr_cell<-thresh(cell_normal, w=10, h=10, offset=0.1)
      colorMode(thr_cell)<-Grayscale
      cmask = opening(thr_cell, kern=makeBrush(7,shape="disc"))
      rm(thr_cell)
      ar<-as.vector(cell_normal)
      ar.sum<-as.numeric(summary(ar))
      open2<-opening(cell_normal>(ar.sum[2]))
      rm(ar)
      colorMode(open2)<-Grayscale
      open_2 = paintObjects(open2,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
      combine<-cmask
      combine[open2 > cmask]<-open2[open2 > cmask]
      combine[gsegg > combine]<-gsegg[gsegg > combine]
      open_2 = paintObjects(combine,toRGB(GFP*150),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
      csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=cmask)
      csegpink <- fillHull(csegpink)
      colorMode(csegpink)<-Grayscale
      csegpink = propagate(cell_normal, gsegg, lambda=1.0e-2, mask=combine)
      csegpink <- fillHull(csegpink)
      colorMode(csegpink)<-Grayscale
      xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
      if (length(xy) < 3 ){
        next_fov()
      }
      if (is.null(nrow(xy))){
        next_fov()
      }
      if (nrow(xy) > 2000 ){
        next_fov()
      }
      cf = computeFeatures.shape(csegpink)
      cf = computeFeatures.shape(csegpink)
      cf = computeFeatures.shape(csegpink)
      cfErea<-data.frame(cf[,1])
      cfErea$num<-row.names(cfErea)
      ci = which(cf[,1] > 35000) 
      csegpink = rmObjects(csegpink, ci,reenumerate = F) 
      rm(cf,ci,cfErea)
      xy.gsegg<-as.numeric(row.names(computeFeatures.moment(gsegg)[,c('m.cx','m.cy')]))
      xy.cseg<- as.numeric(row.names(computeFeatures.moment(csegpink)[,c('m.cx','m.cy')])) 
      ind.deff<-setdiff(xy.gsegg,xy.cseg)
      gsegg<-rmObjects(gsegg,ind.deff,reenumerate=F)
      gsegg<-reenumerate(gsegg)
      csegpink<-reenumerate(csegpink)
      seg_pink = paintObjects(csegpink,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
      display(seg_pink,method="raster")
      display(colorLabels(csegpink),method="raster")
      #Sys.sleep(2)
      xy<-computeFeatures.moment(gsegg)[,c('m.cx','m.cy')]
      if (length(xy) < 3 ){
        next_fov()
      }
      if (is.null(nrow(xy))){
        next_fov()
      }
      if (nrow(xy) > 2000 ){
        next_fov()
      }
      stack<-stackObjects(csegpink,GFP*20,ext = c(200,200))
      px_st = stack
      n = numberOfFrames(px_st)
      px = resize(px_st, 61, 61)
      pxlist = getFrames(px, i=1:n)
      px_A = abind(pxlist, along=0)
      tensor = array_reshape(px_A, c(n, 61, 61, 1))
      preds = predict_on_batch(model, tensor)
      df<-as.data.frame(xy)
      df$Prediction<-as.vector(preds)
      display(GFP*20,"raster")
      celltext = text(x= df[,1], y= df[,2] , labels=round(df[,3],3), col="red", cex = 0.7)
      #Sys.sleep(7)
      ind.na<-which(is.na(df[,3]))
      df[ind.na,3]<-0
      ind<-which(df[,3] < 0.998)
      if (length(ind) < 1 ){
        next_fov()
      }
      #' remove all the cell which were aggragte becouse of missegmentation
      csegpink = rmObjects(csegpink, ind,reenumerate = T)
      gsegg = rmObjects(gsegg, ind,reenumerate = T)
      if (length(table(csegpink)) < 2){
        next_fov()
      }
      stack.pred<-stackObjects(csegpink,GFP*20,ext = c(200,200))
      indd<-which(stack.pred>0)
      stack.pred[indd]<-1
      xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
      if (is.null(nrow(xy))){
        gsegg<-bwlabel(gsegg)
        nseg_bin<-erode(gsegg,makeBrush(17,"line"))
        nseg_bin<-thresh(nseg_bin) ##makes it binery
        display(nseg_bin,method="raster")
        writeImage(x = nseg_bin, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
        closeAllConnections()
        setwd("c:/outproc/")
        while (!sync.txt == 2) { 
          try(
            {
              fileConn<-file("c:/outproc/sync.txt")
              writeLines(c("2"), fileConn)
              close(fileConn)
              print("sync.txt == 2")
              break
            }  
            , silent = TRUE
          )
        } 
        
      }
      df.temp<-as.data.frame(xy)
      df.temp$index<-c(1:dim(stack.pred)[3])
      df.temp$singel_object
      for (m in 1:dim(stack.pred)[3]){
        st.temp<-stack.pred[,,m]
        bg<-length(which(st.temp==0))
        fg<-length(which(st.temp > 0))
        if ((bg/fg) > 10){
          df.temp[m,4]<-"yes"
        } else {
          df.temp[m,4]<-"no"
        }
      } 
      rm(ind)
      ind<-which(df.temp[,4]=="no")
      if (length(ind) < 1 ){
        next_fov()
        
      }
      csegpink = rmObjects(csegpink, ind,reenumerate = T)
      gsegg = rmObjects(gsegg, ind,reenumerate = T)
      seg_call = paintObjects(gsegg,toRGB(GFP*20),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
      display(seg_call,method="raster")
      #Sys.sleep(5)
      gsegg<-bwlabel(gsegg)
      nseg_bin<-erode(gsegg,makeBrush(17,"line"))
      nseg_bin<-thresh(nseg_bin) ##makes it binery
      writeImage(x = nseg_bin, "Binary1.tif", type = "tiff", bits.per.sample = 8L,compression = "LZW" )
      closeAllConnections()
      setwd("c:/outproc/")
      while (!sync.txt == 2) { 
        try(
          {
            
            #closeAllConnections()
            #setwd("c:/outproc/")  
            fileConn<-file("c:/outproc/sync.txt")
            writeLines(c("2"), fileConn)
            close(fileConn)
            print("sync.txt == 2")
            break
          }  
          , silent = TRUE
        )
      }
      
    }  
    , silent = TRUE
  )
  
}
