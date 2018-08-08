source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
biocLite("ggbiplot")
# biocLite("DESeq2")

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")

setwd("/Users/carlykenkel/Dropbox/AIMSpostdoc/PoritesPatch")

 #To upload DESeq package - you need to do this every time you open R
library(DESeq)
library(gplots) # for venn diagram
library(ggplot2)
library(RColorBrewer)
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(plotrix)
library(reshape2)
library(factoextra)


counts=read.table("allcountsPoritesPatchGmapper.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

head(counts) 


length(counts[,1])  #24619 isogroups for Bowtie2; 26578 for Gmapper

names(counts)


readsleft=c()
for (column in names(counts)) {
	val=sum(counts[,column])
	readsleft=append(readsleft,val)}

RLtable=data.frame(cbind(names(counts),readsleft))
RLtable$readsleft=as.numeric(as.character(RLtable$readsleft))
write.csv(RLtable,"readsleft.csv",quote=F) 

#######################Creating table of conditions for your experiment

patch=c(1:length(names(counts)))
patch[grep("B",names(counts))]="Patch"
patch[grep("N",names(counts))]="Normal"
geno=c("X1","X1","X2","X2","X3","X3","X4","X4","X5","X5","X6","X6","X7","X7","X8","X8")


conditions=data.frame(patch,geno)
head(conditions)


############################################original DESeq methods

real=newCountDataSet(counts.nobad1,conditions.nobad1) 
real=estimateSizeFactors(real)

# ####all the data you ever wanted about quality control

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

v="/Users/carlykenkel/Dropbox/AIMSpostdoc/PoritesPatch/DESeq/AQM2"

arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("geno"),force=TRUE) #check .html output file in new folder

###############Remove outliers as detected; repeat arrayQualityMetrics above after regenerating newCountDataSet
###############to confirm that all outliers have been removed 



counts.nobad1=counts[,-c(7,8)] #According to heatmap AQMv1, Sample 4N is serious outlier; 4B questionable. Remove this geno rep and try re-running 
conditions.nobad1=conditions[-c(7,8),]
conditions.nobad1$geno<-factor(conditions.nobad1$geno)

#6B still slightly questionable...but to keep N high so can use gene-est-only, will keep 6B for now

#Are counts low? check
counts2<-counts.nobad1

counts2$low = apply(counts.nobad1[,1:14],1,function(x){sum(x<=10)})
#Sum up samples with counts less than or equal to 10 by isogroup

countsHigh<-counts2[-which(counts2$low>12),] #12 is 85% of 14 samples - how many isogroups have zero counts in 85% or more of samples
nrow(countsHigh)/nrow(counts) #53% of isogroups are sufficiently high 

#NOTE - above is informative, but do NOT pre-filter counts before running through DESeq. The low expression counts are used to calculate the full dispersions model
#Low expression isogroups will be filtered later in the pipeline

###################
### DESEQ 2 #######
###################
library(DESeq2)


dds <- DESeqDataSetFromMatrix(countData = counts.nobad1,
                              colData = conditions.nobad1,
                              design= ~ patch + geno)

dds <- DESeq(dds)

resPatch <- results(dds, contrast=c("patch","Patch","Normal"))

summary(resPatch)

# out of 26553 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 25, 0.094%
# LFC < 0 (down)     : 96, 0.36%
# outliers [1]       : 0, 0%
# low counts [2]     : 18017, 68%
# (mean count < 13)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#meh. I don't love the black box nature of this. would prefer to run the whole analysis myself...

###################
### DESEQ 1 #######
###################

library(DESeq)

real=newCountDataSet(counts.nobad1,conditions.nobad1) 
real=estimateSizeFactors(real) #library size (~total mapped reads)

sizeFactors(real)

real=estimateDispersions(real,method="pooled-CR",sharingMode="gene-est-only",modelFormula=count~geno+patch)  
#can use gene-est-only, IF have >7 reps per treatment combo
# quartz()
plotDispEsts(real)
plotDispEsts(dds) #looks reasonable...but may not be same as regular deseq?

vsd=getVarianceStabilizedData(real)
# write.csv(vsd, file="VSD_allsymgenes_nobadsam_Jan2015.csv", quote=F)

# ######################Determining quality filtering cutoffs

fit0=fitNbinomGLMs(real, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(real, count ~ geno)
fit2=fitNbinomGLMs(real, count ~ patch)

pvals.g<-nbinomGLMTest(fit1,fit0) #testing significance of patch condition
pvals.p<-nbinomGLMTest(fit2,fit0) #testing significance of genotype

pvalue=pvals.g #change to each type of pval: i, t and o
theta=seq(from=0,to=0.8,by=0.02)

filterChoices=data.frame(`mean`=rowMeans(counts(real)),`median`=apply((counts(real)),1,median),`min`=rowMin(counts(real)),`max`=rowMax(counts(real)),`sd`=rowSds(counts(real)))
rejChoices=sapply(filterChoices,function(f) filtered_R(alpha=0.1,filter=f,test=pvalue,theta=theta,method="BH"))
library("RColorBrewer")
myColours=brewer.pal(ncol(filterChoices),"Set1")

# #quartz()
# #windows()
matplot(theta,rejChoices,type="l",lty=1,col=myColours,lwd=2,xlab=expression(theta),ylab="number of rejections")
legend("bottomleft",legend=colnames(filterChoices),fill=myColours)

# #look for peak in graph - corresponds to correct theta and best-fit line for which metric to use - pick best theta for all tests
# #patch=0.5, but squiggly...could make case for 0.2 ; geno=0.2/0.3...lets go with 0.3

# #######################Quality Filtering Data based on theta - get rid of genes with low variance

# #FOR host
rs=rowSds(counts(real)) #using standard deviation as quality filtering metric based on analyses above
theta=0.3 
use=(rs>quantile(rs,probs=theta)) ###
table(use) 
# # use
# use
# FALSE  TRUE 
 # 7978 18600


realFilt=real[use,]
vsd=getVarianceStabilizedData(realFilt)

######################## Now for the real Model Testing

fit0=fitNbinomGLMs(realFilt, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(realFilt, count ~ geno)
fit2=fitNbinomGLMs(realFilt, count ~ patch)
fit3=fitNbinomGLMs(realFilt, count ~ geno+patch)


# testing section
pvals.p<-nbinomGLMTest(fit3,fit1) #testing significance of patch, after accounting for genotype
pvals.g<-nbinomGLMTest(fit1,fit0) #testing significance of genotype


#making non-convergent model p-values NA's - only occur in pvals.p
pvals.p=data.frame(pvals.p)
rownames(fit3)->rownames(pvals.p)
badmods=subset(fit3,(!fit3[,10]))
for ( i in rownames(badmods)){pvals.p[i,1]<-NA}


summary(pvals.p)

pvals.g=data.frame(pvals.g)
rownames(fit1)->rownames(pvals.g)
badmods=subset(fit1,(!fit1[,9]))
for ( i in rownames(badmods)){pvals.g[i,1]<-NA}

summary(pvals.g)

#multiple test correction - adjust p-values using Benjamini-Hochburg

adjp.p<-p.adjust(pvals.p$pvals.p,method="BH")
adjp.g<-p.adjust(pvals.g$pvals.g,method="BH")


do<-(cbind(vsd, "adjp.p" = adjp.p,"pval.p" = pvals.p$pvals.p,"adjp.g" = adjp.g,"pval.g" = pvals.g$pvals.g)) #creating table of all multiple test corrected p-values with variance stabilized count data 


write.csv(do, file="VSDandPVALSnobadsam_deseq1_6aug.csv", quote=F) #writing an output file of vsd plus p-values

########################## counting, venn diagram:

p<-data.frame(do)
patch=row.names(p[p$adjp.p<=0.1 & !is.na(p$adjp.p),])
geno=row.names(p[p$adjp.g<=0.1 & !is.na(p$adjp.g),])


candidates=list("Patch, P<=0.1"=patch,"Geno, P<=0.1"=geno)
library(gplots)
quartz()
venn(candidates)

###################

d2<-read.csv("VSDandPVALSnobadsam_deseq1_6aug.csv")
rownames(d2)<-d2$X

###########Between Groups analysis

library(made4)
 
PatchType=(rep(c("Patch","Normal"),times=7))

tdat<-(d2[,2:15]) 

pca<-bga(tdat,PatchType,type="pca")
plot.bga(pca)
topgenes(pca,n=5)

#NOW, pick out only TOP DEGs
topnum=15 # number of top candidates per comparison. The total number over all comparisons must not exceed the number of samples

GEs<-d2 #rename the dataset

patch=head(GEs[order(GEs$adjp.p),],topnum)
geno=head(GEs[order(GEs$adjp.g),],topnum)
sig=data.frame(rbind(patch[!(patch$X %in% geno$X),],geno)) #remove redundant isogroups
length(sig[,1]) #Top 15 each, none overlapping


#Transpose expression values only
tsig=t(sig[,2:15]) 
tsig[1:5,1:5]


# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE.
library(devtools)

 
pca <- prcomp(tsig,
      center = TRUE,
      scale. = TRUE) 
summary(pca) #First 2 PCs explain 61% of the variation
str(pca)

scores=pca$x

scores[,1:2]

par(mfrow=c(1,2))
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",type="n",main="Patch")
points(scores[c(1,3,5,7,9,11,13),1],scores[c(1,3,5,7,9,11,13),2],pch=1,col="lightgreen")
points(scores[c(2,4,6,8,10,12,14),1],scores[c(2,4,6,8,10,12,14),2],pch=19,col="darkgreen")

plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",type="n",main="Geno")
points(scores[c(1:2),1],scores[c(1:2),2],pch=1,col="pink")
points(scores[c(3:4),1],scores[c(3:4),2],pch=1,col="violet")
points(scores[c(5:6),1],scores[c(5:6),2],pch=1,col="green")
points(scores[c(7:8),1],scores[c(7:8),2],pch=1,col="blue")
points(scores[c(9:10),1],scores[c(9:10),2],pch=1,col="red")
points(scores[c(11:12),1],scores[c(11:12),2],pch=1,col="orange")
points(scores[c(13:14),1],scores[c(13:14),2],pch=1,col="black")

#Picking out top 15 genes in each of the patch and geno categories, shows that patch
#is stronger driver of expression in this instance


#########Generating directional GO output for DESeq results

#patch
fit3$direction=ifelse(fit3$patchPatch>0,1,-1) #red=higher expression in patch
fit3$pval<-(-log((pvals.p$pvals.p+0.0000000001),10))*(fit3$direction)

patch<-cbind("gene"=rownames(fit3),"pval"=fit3$pval) 

write.csv(patch,file="GOpatch.csv",quote=F,row.names=F)

##########################################Categorical GO files for gomwu scripts

d<-d2

d$binary=ifelse(d$adjp.p<=0.1 & !is.na(d$adjp.p),1,0)
d$binary=as.factor(d$binary)
summary(d)
names(d)
GObinary<-d[,c(1,20)]
write.csv(GObinary,file="GObinaryPatch.csv",quote=F,row.names=F)

#And a vsd file for just this "interesting" gene subset

VSDsigOnly<-d[d$binary==1,]
nrow(VSDsigOnly)
head(VSDsigOnly)
write.csv(VSDsigOnly,file="VSDs_GObinaryPatch.csv",quote=F,row.names=F)

####################Now for some heatmaps

gg=read.table("plob_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE) 

d$direction<-fit3$direction

patch=d[d$adjp.p<=0.1 & !is.na(d$adjp.p),] #same genes used in Venn, now to assess direction of change and heatmap plot w/gene names
nrow(patch)
head(patch)

patchup=patch[patch$direction>0 & !is.na(patch$adjp.p),] #genes up-regulated in patch = 83/234
patchdown=patch[patch$direction<0 & !is.na(patch$adjp.p),] #genes down-regulated in patch = 151/234

is=rownames(patchup) #change to list of interest each time
dat=as.data.frame(patchup)

names(dat)
edata=c(2:15) #only columns with your expression data 


#is=c("c_sym_59999", "c_sym_28988", "c_sym_22292", "c_sym_29386" ,"c_sym_71409","c_sym_40040", "c_sym_31320", "c_sym_11928", "c_sym_23149", "c_sym_92871")

#########################append gene names to isogroups
sel=c();gnms=c()
for ( i in is){
	if (i %in% dat$X){
		sel=rbind(sel,dat[dat$X==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}
row.names(sel)=paste(gnms,sel$X,sep=".")

exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp) #should match number of genes in list of interest above

head(exp)
#write.csv(exp,file="patchup_genes.csv",quote=FALSE)


means=apply(exp,1,mean) # means of rows
expc=exp-means #rescale expression data so it's up and down relative to mean

 
library(RColorBrewer)
library(pheatmap)

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=1)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

quartz()
pdf("HeatmapPatchUp.pdf",width=8,height=17)
pheatmap(expc,color=col,cluster_cols=F,clustering_distance_rows="correlation") #plot the heatmap
dev.off()

#oops! also need to reorder heatmap to get patch and normal together
names(expc)
expc<-expc[,c(2,4,6,8,10,12,14,1,3,5,7,9,11,13)] #go up and replot heatmap
