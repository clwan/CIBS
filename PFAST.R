#### PFAST
PFAST_Core<-function(MAT_AR){
  E<-NULL
  S<-order(rowSums(MAT_AR),decreasing = T)
  Cl<-rep(0,nrow(MAT_AR))
  Cr<-rep(0,ncol(MAT_AR))
  Cl[S[1]]<-1
  Cr<-MAT_AR[S[1],]
  for (i in 2:length(S)) {
    Cl1<-Cl
    Cr1<-Cr
    Cl1[S[i]]<-1
    Cr1[which(MAT_AR[S[i],]==0)]<-0
    if(sum(Cr1)<sum(Cr)){
      TEMP<-sum(MAT_AR[S[i],as.logical(Cr1)])-1+sum(Cr-Cr1)-sum(MAT_AR[as.logical(Cl1),as.logical(Cr-Cr1)])
    }else{
      TEMP<-sum(MAT_AR[S[i],])-1
    }
    if(TEMP>0){
      Cl<-Cl1
      Cr<-Cr1
    }else{
      E<-append(E,S[i])
    }
  }
  Core<-list()
  Core[[1]]<-Cl
  Core[[2]]<-Cr
  Core[[3]]<-E
  names(Core)<-c("left","right","Extension")
  return(Core)
}


PFAST_Ext_Core<-function(Core,MAT_A,Thres){
  if(sum(Core[[2]])==1){
    MAT_TEMP<-MAT_A[Core[[3]],as.logical(Core[[2]])]
    Core[[1]][Core[[3]][which(MAT_TEMP>1)]]<-1
  }else{
    MAT_TEMP<-MAT_A[Core[[3]],as.logical(Core[[2]])]
    VEC<-rowSums(MAT_TEMP)
    Core[[1]][Core[[3]][which(VEC>(ncol(MAT_TEMP)*Thres))]]<-1
  }
  
  return(Core)
} 

PFAST<-function(MAT,DIM=1000,Thres){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  
  MAT_A<-MAT
  MAT_B<-NULL
  MAT_C<-NULL
  MAT_AR<-MAT_A
  for (i in 1:DIM) {
    print(sum(MAT_AR))
    if(sum(MAT_AR)<=0.05*sum(MAT)){
      result<-list()
      result[[1]]<-MAT_B
      result[[2]]<-MAT_C
      return(result)
      break()
    }else{
      Core<-PFAST_Core(MAT_AR)
      Core<-PFAST_Ext_Core(Core,MAT_A,Thres)
      MAT_B<-cbind(MAT_B,Core[[1]])
      MAT_C<-rbind(MAT_C,Core[[2]])
      
      TEMP<-Core[[1]] %o% Core[[2]]
      MAT_AR<-MAT_AR-TEMP
      MAT_AR<-apply(MAT_AR>0,2,as.numeric)
    }
  }
  result<-list()
  result[[1]]<-MAT_B
  result[[2]]<-MAT_C
  return(result)
}
