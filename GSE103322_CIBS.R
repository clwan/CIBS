source("PFAST.R")
library(LTMGSCA)
library(Rtsne)
library(mixtools)
library(cluster)
library(RColorBrewer)
library(Seurat)

# Data preparation
load("GSE103322_Cell_key.RData")
File_data<-read.table("GSE103322_TPM.txt",header = T,row.names = 1)
File_data<-File_data[rowSums(File_data)>0,colSums(File_data)>0]


### Expression On and Off by LTMG ###
Zcut_G<-log(Global_Zcut(File_data))
File_LTMG<-LTMG_MAT(MAT = File_data,Zcut_G = Zcut_G,Gene_use = nrow(File_data))
File_state<-File_LTMG[[1]]


#### BMF Analysis ####
File_BM<-File_state>0
BM_use<-PFAST(File_BM,DIM=10,Thres = 0.6)

File_BM<-File_data*(apply(BM_use[[1]][,1:5]%*%BM_use[[2]][1:5,], 2, as.logical))
File_BM<-File_BM[rowSums(File_BM)>0,colSums(File_BM)>0]
File_BM<-log(File_BM+1)


#### Visualization on low dimension. 
TSNE_BM<-Rtsne(t(File_BM),dims = 2, perplexity=30, verbose=TRUE, max_iter = 20000,partial_pca=TRUE)

Cell_BM<-Cell_key[colnames(File_BM)]
Cell_unique<-unique(Cell_key)
cl1<-rainbow(length(Cell_unique))
plot(x=TSNE_BM$Y[which(Cell_BM==Cell_unique[1]),1],
     y=TSNE_BM$Y[which(Cell_BM==Cell_unique[1]),2],
     pch=18,cex=0.7,col=cl1[1],
     xlim=c(min(TSNE_BM$Y)-2,max(TSNE_BM$Y)+2),
     ylim=c(min(TSNE_BM$Y)-2,max(TSNE_BM$Y)+2),,
     frame.plot = F,
     xlab = "t-SNE 1",
     ylab = "t-SNE 2")

