### R code from vignette source 'puma.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: puma.Rnw:67-68
###################################################
options(width = 60)


###################################################
### code chunk number 2: puma.Rnw:114-117 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("puma")


###################################################
### code chunk number 3: puma.Rnw:123-124
###################################################
library(puma)


###################################################
### code chunk number 4: puma.Rnw:128-130 (eval = FALSE)
###################################################
## help(pumaDE)
## ?pumaDE


###################################################
### code chunk number 5: puma.Rnw:134-135 (eval = FALSE)
###################################################
## help(package="puma")


###################################################
### code chunk number 6: puma.Rnw:144-153 (eval = FALSE)
###################################################
## library(oligo)
## estrogenFilenames<-list.celfiles()
## oligo.estrogen<-read.celfiles(celFiles)
## pData(oligo.estrogen) <- data.frame(
##     "estrogen"=c("absent","absent","present","present"
##         ,"absent","absent","present","present")
## ,   "time.h"=c("10","10","10","10","48","48","48","48")
## ,   row.names=rownames(pData(oligo.estrogen))
## )


###################################################
### code chunk number 7: puma.Rnw:158-160
###################################################
library(pumadata)
data(oligo.estrogen)


###################################################
### code chunk number 8: puma.Rnw:162-163
###################################################
show(oligo.estrogen)


###################################################
### code chunk number 9: puma.Rnw:167-168
###################################################
pData(oligo.estrogen)


###################################################
### code chunk number 10: puma.Rnw:177-179 (eval = FALSE)
###################################################
## eset_estrogen_mmgmos<-mmgmos(oligo.estrogen,gsnorm="none")
## eset_estrogen_rma<-rma(oligo.estrogen)


###################################################
### code chunk number 11: puma.Rnw:186-187
###################################################
data(eset_estrogen_mmgmos)


###################################################
### code chunk number 12: puma.Rnw:190-192 (eval = FALSE)
###################################################
## exprs(eset_estrogen_mmgmos)[1,]
## assayDataElement(eset_estrogen_mmgmos,"se.exprs")[1,]


###################################################
### code chunk number 13: puma.Rnw:199-200 (eval = FALSE)
###################################################
## write.reslts(eset_estrogen_mmgmos, file="eset_estrogen_mmgmos.rda")


###################################################
### code chunk number 14: puma.Rnw:207-208 (eval = FALSE)
###################################################
## eset_estrogen_pmmmgmos <- PMmmgmos(oligo.estrogen,gsnorm="none")


###################################################
### code chunk number 15: puma.Rnw:221-225 (eval = FALSE)
###################################################
## setwd(cel.path)   ## go to the directory where you put the CEL files.
## exonFilenames<-c("C0006.CEL","C021.CEL"
##                 ,"F023_NEG.CEL","F043_NEG.CEL") ## or you can use  exonFilenames<-list.celfiles()
## oligo.exon<-read.celfiles(exonFilenames)


###################################################
### code chunk number 16: puma.Rnw:230-231 (eval = FALSE)
###################################################
## eset_gmoExon<-gmoExon(oligo.exon,exontype="Human",gsnorm="none")


###################################################
### code chunk number 17: puma.Rnw:236-237
###################################################
data(eset_gmoExon)


###################################################
### code chunk number 18: puma.Rnw:240-244 (eval = FALSE)
###################################################
## exprs(eset_gmoExon$gene)[1,]
## assayDataElement(eset_gmoExon$gene,"se.exprs")[1,]
## exprs(eset_gmoExon$transcript)[1,]
## assayDataElement(eset_gmoExon$transcript,"se.exprs")[1,]


###################################################
### code chunk number 19: puma.Rnw:251-252 (eval = FALSE)
###################################################
## write.reslts(eset_gmoExon$gene, file="eset_gmoExon_gene")


###################################################
### code chunk number 20: puma.Rnw:261-262 (eval = FALSE)
###################################################
## write.reslts(eset_gmoExon$transcript, file="eset_gmoExon_transcript")


###################################################
### code chunk number 21: puma.Rnw:288-293 (eval = FALSE)
###################################################
## cel.path<-"cel.path"  ## for example  cel.path<-"/home/gao/celData"  
## SampleNameTable<-"SampleNameTable" 
## eset_igmoExon<-igmoExon(cel.path="cel.path",SampleNameTable="SampleNameTable"
##                        , exontype="Human"
##                        , gsnorm="none", condition="Yes")


###################################################
### code chunk number 22: puma.Rnw:304-309 (eval = FALSE)
###################################################
## setwd(cel.path)   ## go to the directory where you put the CEL files.
## library(puma)
## library(oligo)
## oligo.hta<-read.celfiles(celFiles)
## eset_gmhta<-gmhta(oligo.hta,gsnorm="none")


###################################################
### code chunk number 23: puma.Rnw:314-315
###################################################
data(eset_gmhta)


###################################################
### code chunk number 24: puma.Rnw:318-322 (eval = FALSE)
###################################################
## exprs(eset_gmhta$gene)[2,]
## se.exprs(eset_gmhta$gene)[2,]
## exprs(eset_gmhta$transcript)[2,]
## se.exprs(eset_gmhta$transcript)[2,]


###################################################
### code chunk number 25: puma.Rnw:329-330 (eval = FALSE)
###################################################
## write.reslts(eset_gmhta$gene, file="eset_gmhta_gene")


###################################################
### code chunk number 26: puma.Rnw:339-340 (eval = FALSE)
###################################################
## write.reslts(eset_gmhta$transcript, file="eset_gmhta_transcript")


###################################################
### code chunk number 27: puma.Rnw:356-357 (eval = FALSE)
###################################################
## pumapca_estrogen <- pumaPCA(eset_estrogen_mmgmos)


###################################################
### code chunk number 28: puma.Rnw:359-362
###################################################
data(pumapca_estrogen)
data(eset_estrogen_rma)
data(eset_estrogen_mmgmos)


###################################################
### code chunk number 29: puma.Rnw:367-368
###################################################
pca_estrogen <- prcomp(t(exprs(eset_estrogen_rma)))


###################################################
### code chunk number 30: puma.Rnw:374-384
###################################################
par(mfrow=c(1,2))
plot(pumapca_estrogen,legend1pos="right",legend2pos="top",main="pumaPCA")
plot(
	pca_estrogen$x
,	xlab="Component 1"
,	ylab="Component 2"
,	pch=unclass(as.factor(pData(eset_estrogen_rma)[,1]))
,	col=unclass(as.factor(pData(eset_estrogen_rma)[,2]))
,	main="Standard PCA"
)


###################################################
### code chunk number 31: puma.Rnw:397-398 (eval = FALSE)
###################################################
## write.reslts(pumapca_estrogen, file="pumapca_estrogen")


###################################################
### code chunk number 32: puma.Rnw:406-409
###################################################
par(mfrow=c(1,2))
boxplot(data.frame(exprs(eset_estrogen_mmgmos)),main="mmgMOS - No norm")
boxplot(data.frame(exprs(eset_estrogen_rma)),main="Standard RMA")


###################################################
### code chunk number 33: puma.Rnw:423-426
###################################################
eset_estrogen_mmgmos_normd <- pumaNormalize(eset_estrogen_mmgmos)
boxplot(data.frame(exprs(eset_estrogen_mmgmos_normd))
    , main="mmgMOS - median scaling")


###################################################
### code chunk number 34: puma.Rnw:443-444 (eval = FALSE)
###################################################
## eset_estrogen_comb <- pumaComb(eset_estrogen_mmgmos_normd)


###################################################
### code chunk number 35: puma.Rnw:446-447
###################################################
data(eset_estrogen_comb)


###################################################
### code chunk number 36: puma.Rnw:452-453
###################################################
colnames(createContrastMatrix(eset_estrogen_comb))


###################################################
### code chunk number 37: puma.Rnw:462-463 (eval = FALSE)
###################################################
## write.reslts(eset_estrogen_comb, file="eset_estrogen_comb")


###################################################
### code chunk number 38: puma.Rnw:468-470
###################################################
pumaDERes <- pumaDE(eset_estrogen_comb)
limmaRes <- calculateLimma(eset_estrogen_rma)


###################################################
### code chunk number 39: puma.Rnw:475-476 (eval = FALSE)
###################################################
## write.reslts(pumaDERes, file="pumaDERes")


###################################################
### code chunk number 40: puma.Rnw:483-485
###################################################
topLimmaIntGene <- topGenes(limmaRes, contrast=7)
toppumaDEIntGene <- topGenes(pumaDERes, contrast=7)


###################################################
### code chunk number 41: puma.Rnw:493-494
###################################################
plotErrorBars(eset_estrogen_rma, topLimmaIntGene)


###################################################
### code chunk number 42: puma.Rnw:509-510
###################################################
plotErrorBars(eset_estrogen_mmgmos_normd, topLimmaIntGene)


###################################################
### code chunk number 43: puma.Rnw:523-524
###################################################
plotErrorBars(eset_estrogen_mmgmos_normd, toppumaDEIntGene)


###################################################
### code chunk number 44: puma.Rnw:547-550
###################################################
data(eset_mmgmos)
eset_mmgmos_100 <- eset_mmgmos[1:100,]
pumaCombImproved <- pumaCombImproved(eset_mmgmos_100)


###################################################
### code chunk number 45: puma.Rnw:555-556
###################################################
colnames(createContrastMatrix(pumaCombImproved))


###################################################
### code chunk number 46: puma.Rnw:563-564
###################################################
write.reslts(pumaCombImproved,file="eset_mmgmo_combimproved")


###################################################
### code chunk number 47: puma.Rnw:569-571
###################################################

pumaDEResImproved <- pumaDE(pumaCombImproved)


###################################################
### code chunk number 48: puma.Rnw:576-577
###################################################
write.reslts(pumaDEResImproved,file="pumaDEResImproved")


###################################################
### code chunk number 49: puma.Rnw:587-588
###################################################
pumaClust_estrogen <- pumaClust(eset_estrogen_mmgmos, clusters=7)


###################################################
### code chunk number 50: puma.Rnw:593-595
###################################################
summary(as.factor(pumaClust_estrogen$cluster))
pumaClust_estrogen$centers


###################################################
### code chunk number 51: puma.Rnw:597-598 (eval = FALSE)
###################################################
## write.csv(pumaClust_estrogen$cluster, file="pumaClust_clusters.csv")


###################################################
### code chunk number 52: puma.Rnw:609-617
###################################################
data(Clustii.exampleE)
data(Clustii.exampleStd)
  r<-vector(mode="integer",0)
  for (i in c(1:20))
    for (j in c(1:4))
      r<-c(r,i)
cl<-pumaClustii(Clustii.exampleE,Clustii.exampleStd,
                mincls=6,maxcls=6,conds=20,reps=r,eps=1e-3)


###################################################
### code chunk number 53: puma.Rnw:635-636 (eval = FALSE)
###################################################
## oligo.estrogen<-read.celfiles(celFiles,pkgname="pd.hg.u95a")


###################################################
### code chunk number 54: puma.Rnw:643-644
###################################################
table(table(oligo::probeNames(oligo.estrogen)))


###################################################
### code chunk number 55: puma.Rnw:677-680 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("snow")


###################################################
### code chunk number 56: puma.Rnw:687-689 (eval = FALSE)
###################################################
## library(snow)
## cl <- makeCluster(c("node01", "node02", "node03", "node04"), type = "SOCK")


###################################################
### code chunk number 57: puma.Rnw:694-703 (eval = FALSE)
###################################################
## library(puma)
## data(affybatch.estrogen)
## pData(affybatch.estrogen) <- data.frame(
## 	"level"=c("twenty","twenty","ten")
## ,	"batch"=c("A","B","A")
## ,	row.names=rownames(pData(affybatch.example)))
## eset_mmgmos <- mmgmos(oligo.estrogen)
## system.time(eset_comb_1 <- pumaComb(eset_mmgmos))
## system.time(eset_comb_4 <- pumaComb(eset_mmgmos, cl=cl))


###################################################
### code chunk number 58: puma.Rnw:707-709 (eval = FALSE)
###################################################
## library(snow)
## cl <- makeCluster(c("localhost", "localhost"), type = "SOCK")


###################################################
### code chunk number 59: puma.Rnw:716-719 (eval = FALSE)
###################################################
## library(snow)
## cl <- makeCluster(c("node01", "node01", "node02", "node02"
## , "node03", "node03", "node04", "node04"), type = "SOCK")


###################################################
### code chunk number 60: puma.Rnw:726-728 (eval = FALSE)
###################################################
## library(snow)
## cl<-makeCluster(c("node01","node02","node03","node04"),type = "SOCK")


###################################################
### code chunk number 61: puma.Rnw:733-737 (eval = FALSE)
###################################################
## library(puma)
## library(oligo)
## object<-read.celfiles("filename.CEL")
## eset<-gmoExon(object,exontype="Human",GT="gene",gsnorm="none",cl=cl)


###################################################
### code chunk number 62: puma.Rnw:742-744 (eval = FALSE)
###################################################
## ibrary(snow)
## cl<-makeCluster(c("loaclhost","localhost"),type = "SOCK")


###################################################
### code chunk number 63: puma.Rnw:749-752 (eval = FALSE)
###################################################
## library(snow)
## cl<-makeCluster(c("node01","node01","node02" ,"node02"
## ,"node03" ,"node03" ,"node04" ,"node04"), type = "SOCK")


###################################################
### code chunk number 64: puma.Rnw:759-761 (eval = FALSE)
###################################################
## library(snow)
## cl<-makeCluster(c("node01","node02","node03","node04"),type = "SOCK")


###################################################
### code chunk number 65: puma.Rnw:766-770 (eval = FALSE)
###################################################
## library(puma)
## library(oligo)
## object<-read.celfiles("filename.CEL")
## eset<-gmhta(object,gsnorm="none",cl=cl)


###################################################
### code chunk number 66: puma.Rnw:775-777 (eval = FALSE)
###################################################
## ibrary(snow)
## cl<-makeCluster(c("loaclhost","localhost"),type = "SOCK")


###################################################
### code chunk number 67: puma.Rnw:782-785 (eval = FALSE)
###################################################
## library(snow)
## cl<-makeCluster(c("node01","node01","node02" ,"node02"
## ,"node03" ,"node03" ,"node04" ,"node04"), type = "SOCK")


###################################################
### code chunk number 68: puma.Rnw:806-810 (eval = FALSE)
###################################################
## library(Rmpi)
## library(snow)
## cl <- makeCluster(4)
## clusterEvalQ(cl, library(puma))


###################################################
### code chunk number 69: puma.Rnw:820-821
###################################################
sessionInfo()


