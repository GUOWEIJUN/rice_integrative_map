#####
setwd("./")
## ----message=FALSE------------------------------------------------------------
library(methylKit)
rm(list=ls())
############
#####CG
file.list=list("../../data/CG/CK-ES.CG-L_R1.txt","../../data/CG/CK-ES.CG-L_R2.txt","../../data/CG/DR-ES.CG-L_R1.txt ","../../data/CG/DR-ES.CG-L_R2.txt ")

###直接进行差异甲基化分析
myobj=methRead(file.list,
               sample.id=list("CK1","CK2","Treat1","Treat2"),
               assembly="rice",
               treatment=c(0,0,1,1),
               context="CpG",
               mincov = 4
)            #2022年NP选取的是4倍

# myobj=methRead(file.list,
#                sample.id=list("CK1","CK2","Treat1","Treat2","Treat3","Treat4",
#                               "Treat5","Treat6","Treat7","Treat8","Treat9","Treat10"),
#                assembly="rice",
#                treatment=c(0,0,1,1,1,1,1,1,1,1,1,1),
#                context="CpG",
#                mincov = 4
# )            #2022年NP选取的是4倍

tiles = tileMethylCounts(myobj,win.size=200,step.size=200,cov.bases = 10)
meth_tiles = unite(tiles)
myDiff_tiles = calculateDiffMeth(meth_tiles)  #,mc.cores=8 'mc.cores' > 1 is not supported on Windows

write.table(myDiff_tiles,file = "all_info.txt",sep = "\t",quote = F)

#提取差异区域
myDiff10p.hyper=getMethylDiff(myDiff_tiles,difference=10,qvalue=0.05,type="hyper")
myDiff10p.hypo=getMethylDiff(myDiff_tiles,difference=10,qvalue=0.05,type="hypo")
myDiff10p=getMethylDiff(myDiff_tiles,difference=10,qvalue=0.05)

write.table(myDiff10p.hyper,file = "hyper.CG.txt",sep = "\t",quote = F)
write.table(myDiff10p.hypo,file = "hypo.CG.txt",sep = "\t",quote = F)
write.table(myDiff10p,file = "alldiff.txt",sep = "\t",quote = F)

##注释文件
library(ChIPseeker)   #加载R包
library(RColorBrewer)
library(GenomicFeatures)   #这个包主要有一个函数将gff文件转成背景注释文件。
rm(list=ls())

NIP_db=makeTxDbFromGFF("../NIP.gtf",format="gtf",taxonomyId=4530)

myDiff10p.hyper <- readPeakFile("hyper.CG.txt",as="GRanges")
myDiff10p.hypo <- readPeakFile("hypo.CG.txt",as="GRanges")



#----1.Peak Annotation------------------------------------------------------

hyper <- annotatePeak(myDiff10p.hyper, tssRegion=c(-2000, 2000),TxDb=NIP_db)
hypo <-  annotatePeak(myDiff10p.hypo, tssRegion=c(-2000, 2000),TxDb=NIP_db)

#---------2.Visualize Genomic Annotation------------------------------------
#2.2 将每个peak注释的基因信息写入文件
write.table(hyper@anno, file="myDiff10p.hyper_annotation.txt", sep="\t", quote=F, row.names=F,col.names = T)
write.table(hypo@anno, file="myDiff10p.hypo_annotation.txt", sep="\t", quote=F, row.names=F,col.names = T)








































# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("CK1","CK2","Treat1","Treat2"),
           assembly="rice",
           treatment=c(1,1,0,0),
           context="CpG",
           mincov = 10
           )

## -----------------------------------------------------------------------------
#getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

## -----------------------------------------------------------------------------
#绘图，甲基化数据统计，可以对四个对象都进行绘图
png("CK-ES_CG-1.png",width = 800,height = 600)
getMethylationStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
dev.off()
png("CK-ES_CG-2.png",width = 800,height = 600)
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
dev.off()
png("DR-ES_CG-1.png",width = 800,height = 600)
getMethylationStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
dev.off()
png("DR-ES_CG-2.png",width = 800,height = 600)
getMethylationStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
dev.off()

## -----------------------------------------------------------------------------
#覆盖度
png("CK-ES_CG.coverage-1.png",width = 800,height = 600)
getCoverageStats(myobj[[1]],plot=TRUE,both.strands=FALSE)
dev.off()
png("CK-ES_CG.coverage-2.png",width = 800,height = 600)
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
dev.off()
png("DR-ES_CG.coverage-1.png",width = 800,height = 600)
getCoverageStats(myobj[[3]],plot=TRUE,both.strands=FALSE)
dev.off()
png("DR-ES_CG.coverage-2.png",width = 800,height = 600)
getCoverageStats(myobj[[4]],plot=TRUE,both.strands=FALSE)
dev.off()



## -----------------------------------------------------------------------------
# filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
#                                       hi.count=NULL,hi.perc=99.9)
# lo.count:指定最小深度
# hi.count:指定最大深度
# lo.perc:
# hi.perc:
## -----------------------------------------------------------------------------
#合并样本数据：
meth=unite(myobj, destrand=FALSE)

## -----------------------------------------------------------------------------
head(meth)

## ----eval=FALSE---------------------------------------------------------------
#  # creates a methylBase object,
#  # where only CpGs covered with at least 1 sample per group will be returned
#  
#  # there were two groups defined by the treatment vector,
#  # given during the creation of myobj: treatment=c(1,1,0,0)
#  meth.min=unite(myobj,min.per.group=1L)

## -----------------------------------------------------------------------------
#样本的相关性：getCorrelation
png("CK-correlation.png",width = 800,height = 600)
getCorrelation(meth,plot=TRUE)
dev.off()
## -----------------------------------------------------------------------------
#样本聚类：
png("DR-ES-cluster.png",width = 800,height = 600)
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
dev.off()
## ----message=FALSE------------------------------------------------------------
png("DR-ES-cluster-hc.png",width = 800,height = 600)
hc = clusterSamples(meth, dist="correlation", method="ward.D2", plot=FALSE)
dev.off()
## -----------------------------------------------------------------------------
png("DR-ES-PCA-barplot.png",width = 800,height = 600)
PCASamples(meth, screeplot=TRUE)
dev.off()
## -----------------------------------------------------------------------------
png("DR-ES-PCA.png",width = 800,height = 600)
PCASamples(meth)
dev.off()
## -----------------------------------------------------------------------------
# make some batch data frame
# this is a bogus data frame
# we don't have batch information
# for the example data
#批量效应用于检查哪些主要成分与潜在的批处理效果上统计相关。
#如果某些主要成分是批次效应造成的结果，就需要使用removeComp将其删除。
# sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
#                             age=c(19,34,23,40))
# 
# as=assocComp(mBase=meth,sampleAnnotation)
# as
# 
# construct a new object by removing the first pricipal component
# from percent methylation value matrix
# newObj=removeComp(meth,comp=1)

## -----------------------------------------------------------------------------
# mat=percMethylation(meth)
# 
# # do some changes in the matrix
# # this is just a toy example
# # ideally you want to correct the matrix
# # for batch effects
# mat[mat==100]=80
#  
# # reconstruct the methylBase from the corrected matrix
# newobj=reconstruct(mat,meth)
# 
# ## ----warning=FALSE------------------------------------------------------------



# 3.1 单碱基差异分析
# 在methylKit的差异分析针对的是合并后的甲基化表达谱，
# 上面的合并文件中每一行是一个甲基化位点，那么差异分析的结果就是差异甲基化位点，
# 使用函数calculateDiffMeth执行差异分析。

## -----------------------------------------------------------------------------
#myDiff=calculateDiffMeth(meth,mc.cores=8)   #多核并行操作

## -----------------------------------------------------------------------------
# myDiff25p = getMethylDiff(myDiff,difference = 25,qvalue = 0.01,type = "hyper" / "hypo" / "all")
# 
# hyper表示相比control组，treatment组中的甲基化C更多；
# hypo则相反，表示treatment组中的甲基化C比control组中少。


# get hyper methylated bases
# myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")  #超甲基化
# #
# write.table(myDiff25p.hyper,file = "DR-ESvsCK_ES.hyper.CG.txt",sep = "\t",quote = F)
# # get hypo methylated bases
# myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
# 
# write.table(myDiff25p.hypo,file = "DR-ESvsCK_ES.hypo.CG.txt",sep = "\t",quote = F)
# 
# 
# #
# #
# # get all differentially methylated bases
# myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
# write.table(myDiff25p.hypo[,-1],file = "DR-ESvsCK_ES.alldiff.CG.txt",sep = "\t",quote = F)

## -----------------------------------------------------------------------------
# diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=20)

## ----eval=FALSE---------------------------------------------------------------
#  
#  sim.methylBase1<-dataSim(replicates=6,sites=1000,
#                           treatment=c(rep(1,3),rep(0,3)),
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  
#  my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
#                                  overdispersion="MN",test="Chisq",mc.cores=1)


# 3.2 区域差异分析
# 执行区域分析之前，首先要合并甲基化区域文件，区域比较时，对单个碱基的覆盖度要求较低。
# 设定窗口和步长，比较区域内甲基化差异
# myobj_lowCov = methRead(file.list,
#            sample.id=list("test1","test2","ctrl1","ctrl2"),
#            assembly="hg18",
#            treatment=c(1,1,0,0),
#            context="CpG",
#            mincov = 3
#            )
# myobj_new=methRead(file.list,
#                sample.id=list("H8_ES_1","H8_ES_2","CK_ES_1","CK_ES_2"),
#                assembly="rice",
#                treatment=c(1,1,0,0),
#                context="CpG",
#                mincov = 10
# )
# 平铺窗口分析：tileMethylCounts
# 差异甲基化区域的分析:按照滑动窗口的方式定义甲基化区域，默认窗口大小为10000 bp ，步长为10000bp.
# 3.2 区域差异分析
# 执行区域分析之前，首先要合并甲基化区域文件，区域比较时，对单个碱基的覆盖度要求较低。
# 设定窗口和步长，比较区域内甲基化差异
myobj=methRead(file.list,
               sample.id=list("CK_ES_CG_1","CK_ES_CG_2","DR-ES_CG_1","DR-ES_CG_2"),
               assembly="rice",
               treatment=c(1,1,0,0),
               context="CpG",
               mincov = 4
)            #2022年NP选取的是4倍

tiles = tileMethylCounts(myobj,win.size=200,step.size=200,cov.bases = 10)
#head(tiles[[1]],3)
meth_tiles = unite(tiles)
#head(meth_tiles)

#读取文件
#meth_tiles = readMethylDB(mydbpath)
#重新选择样品
#myobjDB2 = reorganize(myobjDB,sample.ids=c("test1","ctrl2"),treatment=c(1,0),suffix = "output_name")

myDiff_tiles = calculateDiffMeth(meth_tiles)  #,mc.cores=8 'mc.cores' > 1 is not supported on Windows
#提取差异区域
myDiff20p.hyper=getMethylDiff(myDiff_tiles,difference=20,qvalue=0.05,type="hyper")
myDiff20p.hypo=getMethylDiff(myDiff_tiles,difference=20,qvalue=0.05,type="hypo")
myDiff20p=getMethylDiff(myDiff_tiles,difference=20,qvalue=0.05)


write.table(myDiff25p.hyper,file = "DMR_DR-ESvsCK_ES.hyper.CG.txt",sep = "\t",quote = F)

write.table(myDiff25p.hypo,file = "DMR_DR-ESvsCK_ES.hypo.CG.txt",sep = "\t",quote = F)

write.table(myDiff20p,file = "DMR_DR-ESvsCK_ES.alldiff.CG.txt",sep = "\t",quote = F)

diffMethPerChr(myDiff20p,plot=TRUE,qvalue.cutoff=0.05, meth.cutoff=20)

## ----eval=FALSE---------------------------------------------------------------
#  
#  covariates=data.frame(age=c(30,80,34,30,80,40))
#  sim.methylBase<-dataSim(replicates=6,sites=1000,
#                          treatment=c(rep(1,3),rep(0,3)),
#                          covariates=covariates,
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  my.diffMeth3<-calculateDiffMeth(sim.methylBase,
#                                  covariates=covariates,
#                                  overdispersion="MN",test="Chisq",mc.cores=1)

## ---- eval=FALSE--------------------------------------------------------------
 myDiff=calculateDiffMeth(meth,mc.cores=2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")

library(ChIPseeker)   #加载R包
library(RColorBrewer)
library(GenomicFeatures)   #这个包主要有一个函数将gff文件转成背景注释文件。
rm(list=ls())

#covplot(peak, weightCol=5)  #这个图是Peaks在染色体上的分布。



NIP_db=makeTxDbFromGFF("../../NIP.gtf",format="gtf",taxonomyId=4530)

myDiff25p.hyper <- readPeakFile("DR-ESvsCK_ES.hyper.CG.txt",as="GRanges")
myDiff25p.hypo <- readPeakFile("DR-ESvsCK_ES.hypo.CG.txt",as="GRanges")



#----1.Peak Annotation------------------------------------------------------

hyper <- annotatePeak(myDiff25p.hyper, tssRegion=c(-2000, 2000),TxDb=NIP_db)
hypo <-  annotatePeak(myDiff25p.hypo, tssRegion=c(-2000, 2000),TxDb=NIP_db)

#---------2.Visualize Genomic Annotation------------------------------------


#2.2 将每个peak注释的基因信息写入文件
write.table(hyper@anno, file="myDiff25p.hyper_annotation.txt", sep="\t", quote=F, row.names=F,col.names = T)
write.table(hypo@anno, file="myDiff25p.hypo_annotation.txt", sep="\t", quote=F, row.names=F,col.names = T)



#2.3 可视化
#饼图
png('pie.png',w=2000,h=1000,res=300,units="px")
plotAnnoPie(peakAnno_R1_466@anno)
dev.off()
#Bar Plot
png('Bar.png',w=2000,h=1000,res=300,units="px")
plotAnnoBar(peakAnno_R1_466@anno)
dev.off()
#重叠图
png('vennpie.png',w=2000,h=1000,res=300,units="px")
vennpie(peakAnno_R1_466@anno)
dev.off()
#重叠图
png('upsetplot.png',w=2000,h=1000,res=300,units="px")
upsetplot(peakAnno)
dev.off()

png('upsetplot1.png',w=2000,h=1000,res=300,units="px")
upsetplot(peakAnno,vennpie=TRUE)
dev.off()
#2.4
promoter <- getPromoters(TxDb=NIP_db, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
#查看数据维数
dim(tagMatrix)
#[1] 34131  6001
#将每列求和
a <- colSums(tagMatrix)
#输出文件
write.table(a, file="../H3K4me2.txt", sep="\t", quote=F, row.names=F,col.names = F)

#2.5 Visualize distribution of TF-binding loci relative to TSS-------
png('Distribution of transcription factor-binding loci relative to TSS.png',w=2000,h=1000,res=300,units="px")
plotDistToTSS(peakAnno)
dev.off()

#2.6 Heatmap of ChIP binding to TSS regions
png('Heatmap of ChIP binding to TSS regions.png',w=1000,h=5000,res=300,units="px")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()
#（或者peakHeatmap(suz12, TxDb=txdb, upstream=3000, downstream=3000, color="blue")）
#2.7 Average Profile of ChIP peaks binding to TSS region(Running bootstrapping for tag matrix...		 2020-03-07 run 16:32:54 )
png('Average Profile of ChIP peaks binding to TSS region.png',w=2000,h=1000,res=300,units="px")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf=0.95,resample = 1000,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()
#2.8 Functional enrichment analysis-----------------------
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)

dotplot(pathway2)


  
