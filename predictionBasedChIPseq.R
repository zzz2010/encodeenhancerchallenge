.libPaths("/home/ubuntu/R/x86_64-pc-linux-gnu-library/3.1/")

getRegionW<-function(r){
	ss=unlist(strsplit(r,",",fixed=T))
	as.integer(ss[length(ss)])-as.integer(ss[length(ss)-1])
}
getRegionW2<-function(r){
        ss=unlist(strsplit(r,"|",fixed=T))
        as.integer(ss[length(ss)])-as.integer(ss[length(ss)-1])
}

FeatureMat=read.table("data/FeatureMat.txt",header=T,sep="\t")
ResponseMat=read.table("data/ResponseMat.txt",header=T,sep="\t")

nFeat=ncol(FeatureMat)

regionWidths.train=unlist(lapply(rownames(ResponseMat),getRegionW))
testFeatureMat=read.table("data/LBNLtestFeatureMat.txt",header=T,sep="\t")
regionWidths.test=unlist(lapply(rownames(testFeatureMat),getRegionW2))
combineFeatureMat=rbind(FeatureMat,testFeatureMat)
weights=ResponseMat
weights[weights<0.1]=-log10(weights[weights<0.1])
print(weights[1:100,1])
dim(FeatureMat)
dim(ResponseMat)
dim(testFeatureMat)
length(regionWidths.test)
##impute missing value##
combineFeatureMat[combineFeatureMat==-1]=NA
cacheFn="/tmp/imputed_result2.RData"
if(!file.exists(cacheFn)){
library(mi)
mdf=missing_data.frame(combineFeatureMat)
show(mdf)
nchain=24
options(mc.cores = nchain)
imputations <- mi(mdf, n.iter = 30, n.chains = nchain, max.minutes = 20)
round(mipply(imputations, mean, to.matrix = TRUE), 3)
dfs <- complete(imputations, m = nchain)
save(dfs,file=cacheFn)
}else{
load(cacheFn)
}
imputedFeatureMat.combine=Reduce("+",dfs)/length(dfs)
imputedFeatureMat=imputedFeatureMat.combine[1:nrow(FeatureMat),]
imputedFeatureMat.test=imputedFeatureMat.combine[(nrow(FeatureMat)+1):nrow(imputedFeatureMat.combine),]
###confirm mi is better impute.knn##
##library(impute)
#imputedFeatureMat <- impute.knn(as.matrix(FeatureMat))$data
#############

library(glmnet)
require(doMC)
nfolds=24
registerDoMC(cores=nfolds)

getProbs<-function(k){
fit=cv.glmnet(as.matrix(cbind(imputedFeatureMat[,1:nFeat],regionWidths.train)),ResponseMat[,k]>0.1,alpha=0, parallel = TRUE,family = "binomial",nfolds=nfolds,type.measure="auc" ,weights=weights[,k])
bestI=which.max(fit$cvm)

print(paste(colnames(ResponseMat)[k]," AUC:",fit$cvm[bestI]))
testMat=cbind(imputedFeatureMat.test[,1:nFeat],regionWidths.test)
predict(fit$glmnet.fit,as.matrix(testMat),s=fit$lambda[bestI],type='response')

}
brain=getProbs(1)
heart=getProbs(2)
other=getProbs(3)
anyE=1-(1-brain)*(1-heart)*(1-other)
prediction=cbind(anyE,brain,heart)
colnames(prediction)<-c("any","brain","heart")
rownames(prediction)<-rownames(testFeatureMat)

write.table(prediction,file="PredictionUsingChIPSeq.txt",sep="\t",quote=F)


