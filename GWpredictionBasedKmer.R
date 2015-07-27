.libPaths("/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.1/")
library(Matrix)
labels=read.table("data/vistaLabels.txt",header=F)
Ymat=matrix(0,max(as.numeric(labels[,1])),3)

seli=as.integer(labels[labels[,2]=="forebrain",1])
Ymat[seli,1]=1

seli=as.integer(labels[labels[,2]=="heart",1])
Ymat[seli,2]=1

seli=as.integer(labels[labels[,2]=="any",1])
Ymat[seli,3]=1
colnames(Ymat)=c("forebrain","heart","any")
fgMat=read.table("data/vista.fa.smat_8mer.gz",header=F)

fgMat2=read.table("data/vista.fa.5merPair.smat.gz",header=F)
trainMat=sparseMatrix(i=c(fgMat[,1],fgMat2[,1]),j=c(fgMat[,2],fgMat2[,2]+max(fgMat[,2])+1),x=c(fgMat[,3],fgMat2[,3]),index1=F)

###make prediction
testfgMat=read.table("data/GW5k.mm10.fa.smat_8mer.gz",header=F)
testfgMat2=read.table("data/GW5k.mm10.fa.5merPair.vistasel.smat.gz",header=F)
testMat=sparseMatrix(i=c(testfgMat[,1],testfgMat2[,1]),j=c(testfgMat[,2],testfgMat2[,2]+max(testfgMat[,2])+1),x=c(testfgMat[,3],testfgMat2[,3]),index1=F)


selFeat=1:ncol(trainMat)
library(glmnet)
library(class)
library(randomForest)
library(pROC)
require(doMC)
nfolds=5
registerDoMC(cores=nfolds)
getProbs<-function(k){
cacheFn=paste("data/kmermodel",k,".RData",sep="")
if(!file.exists(cacheFn)){
filter=(rowSums(Ymat)==0|Ymat[,k]==1)
fit=cv.glmnet(trainMat[filter,selFeat],Ymat[filter,k],nfolds=nfolds,parallel=T,alpha=0,family="binomial",type.measure="auc")
bestI=which.max(fit$cvm)
print(paste(colnames(Ymat)[k] ,fit$cvm[bestI]))
save(fit,file=cacheFn)
}else{
load(cacheFn)
}
bestI=which.max(fit$cvm)
predict(fit$glmnet.fit,testMat,s=fit$lambda[bestI],type='response')
}



brain=getProbs(1)
heart=getProbs(2)
anyE=getProbs(3)
prediction=cbind(anyE,brain,heart)
colnames(prediction)<-c("any","brain","heart")
rownames(prediction)<-read.table("data/GW5k.mm10.fa.name",stringsAsFactors=F)[,1]
o=order(-prediction[,1])
write.table(prediction[o,],file="GenomeWide5kPredictionUsingKmer.txt",sep="\t",quote=F)


