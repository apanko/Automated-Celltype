
JoinedZhang<-read.csv("Zhang_FullDataset.csv", row.names=1, header=T)


#Entering in cell identities:
colnames(JoinedZhang_AsNumMatrix_Log2)
ExperimentalCellType<-c("Astrocyte", "Astrocyte", "Neuron", "Neuron", "OPC", "OPC", "NFO", "NFO", "Oligodendrocyte", "Oligodendrocyte", "Microglia", "Microglia", "Endothelial", "Endothelial", "WholeBrain", "WholeBrain", "WholeBrain")

table(ExperimentalCellType)

cbind(ExperimentalCellType, colnames(JoinedZhang_AsNumMatrix_Log2))



#Taking a look at the Zhang data:
is.numeric(JoinedZhang)
is.numeric(as.matrix(JoinedZhang))
#That works.

GeneNamesForJoinedZhang<-row.names(JoinedZhang)

JoinedZhang_AsNumMatrix<-as.matrix(JoinedZhang)
is.numeric(JoinedZhang_AsNumMatrix)

#Let's take a peek at their distribution:

#I wonder if the principal components of variation in the data correlate with those variables?  

#I should probably log transform the reads first before messing around with PCA:
#Let's check out the average distribution of the signal first:
temp<-apply(JoinedZhang_AsNumMatrix, 1, mean)
png("Histogram_AverageReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
max(temp)
#14317.77
median(temp)
[1] 1.996951

temp<-apply(JoinedZhang_AsNumMatrix, 1, max)
png("Histogram_MaxReadsPerProbe.png")
hist(temp, breaks=1000)
dev.off()
#Wow, that is super skewed with a huge number of probes with an average, max, or median of 0. Yeah, let's log transform that data.

#Oh wait - log transformation in RNAseq data is awkward, because there are 0's, which convert into -Inf. It seems that folks typically log transform shifted data instead (data+1)
JoinedZhang_AsNumMatrix_Log2<-log2((JoinedZhang_AsNumMatrix+1))
row.names(JoinedZhang_AsNumMatrix_Log2)<-GeneNamesForJoinedZhang
write.csv(JoinedZhang_AsNumMatrix_Log2, "JoinedZhang_AsNumMatrix_Log2.csv")

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, mean)
png("Histogram_AverageReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, max)
png("Histogram_MaxReadsPerProbeLog2.png")
hist(temp, breaks=1000)
dev.off()

#That's still super skewed, but not as bad as before. I'm guessing that it is difficult to get reads >0 for many transcripts in single cell type data. That means that any analysis based on sample-sample correlations is going to be largely driven by a few highly-expressed transcripts.
#Apparently there are other transformation methods out there that work better for RNAseq - VST is implemented by DESeq. VOOM applies a tranformation to the read counts that supposedly then makes the RNAseq data compatible with analyses included in the limma package. https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/


#Time to recycle some code:
library(gdata)
library(fields)
library(stats)
library(car)
library(affy)
library(preprocessCore)
library(multtest)

#The fact that this data (unlike our microarray data) isn't quantile normalized may make the results of PCA a little funky. Let's start with sample-sample correlations:

#Visualize the sample-sample correlations using a heatmap:

temp<-cor(JoinedZhang_AsNumMatrix_Log2)
colnames(temp)<-colnames(JoinedZhang_AsNumMatrix_Log2)
row.names(temp)<-colnames(JoinedZhang_AsNumMatrix_Log2)

png("09 Sample Sample Correlations Heatmap.png")
heatmap(temp, main="Visualizing correlations between entire samples", xlab="Red=Less correlated, Light yellow=Highly correlated")
dev.off()
#Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)

#Visualize the sample-sample correlations using a boxplot:
png("09 Boxplot Sample Sample Correlations.png", width=2000, height=300)
boxplot(data.frame(cor(JoinedZhang_AsNumMatrix_Log2)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
Median10thQuantile<-median(apply((cor(JoinedZhang_AsNumMatrix_Log2)), 1, quantile, 0.1))
MedianQuantile<-median(apply((cor(JoinedZhang_AsNumMatrix_Log2)), 1, quantile, 0.5))
abline(a=Median10thQuantile, b=0, col=2)
abline(a=MedianQuantile, b=0, col=3)
mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
dev.off()

#It might be worthwile to perhaps try running things with all the genes with extremely low reads thrown out:

temp<-apply(JoinedZhang_AsNumMatrix_Log2, 1, max)
sum(temp<2)
[1] 9787
#That's a lot of genes with 2 reads (log2=1) or less to begin with.  
sum(temp<3)
[1] 11889
#And almost half of them don't meet a standard cut off of 8 reads. Ouch.
sum(temp<3)/length(temp)
#[1] 0.5292939
#...or 53%. Dang. That's worse than the single-cell data. Does that mean that the read depth was shallow?



#Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
pcaNormFiltered<-prcomp(t(JoinedZhang_AsNumMatrix_Log2))
tmp<-pcaNormFiltered$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")

PC1<-pcaNormFiltered$x[,1]
PC2<-pcaNormFiltered$x[,2]

PC3<-pcaNormFiltered$x[,3]
PC4<-pcaNormFiltered$x[,4]

#Output a scree plot for the PCA:
png("09 PCA Scree Plot1.png")
plot(summary(pcaNormFiltered)$importance[2,]~(c(1:ncol(JoinedZhang_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("09 PCA Scree Plot2.png")
plot(summary(pcaNormFiltered)$importance[3,]~(c(1:ncol(JoinedZhang_AsNumMatrix_Log2))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("09 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(ExperimentalCellType))
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("09 PC3 vs PC4.png")
plot(PC3~PC4, main="Principal Components Analysis of Normalized Filtered Data", col=as.factor(ExperimentalCellType))
dev.off()

#Nice clustering.

###############################################################
#I'm going to z-score my log-transformed data.
#I wonder where the NAs are coming from.
#Ah -after messing around, I believe either NaNs or Inf are produced whenever a gene doesn't have any variability associated with it
#So let's remove the completely invariable data first:

JoinedZhang_StDev<-apply(JoinedZhang_AsNumMatrix_Log2, 1, sd) 
sum(JoinedZhang_StDev==0)
#[1] 5314
#Dang....

JoinedZhang_AsNumMatrix_Log2_NoSD0<-JoinedZhang_AsNumMatrix_Log2[JoinedZhang_StDev>0,]
temp<-GeneNamesForJoinedZhang[-c(22086:22088)]
GeneNamesForJoinedZhang_NoSD0<-temp[JoinedZhang_StDev>0]


ZscoreZhang<-t(scale(t(JoinedZhang_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
write.csv(ZscoreZhang, "ZscoreZhang.csv")
sum(is.na(ZscoreZhang))
#looks good now

###############################################################

#My newer version of the cell type analysis: stolen from the ABA code

CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)

colnames(CellTypeSpecificGenes_Master3)

[1] "Umbrella.Cell.Type"    "Specific.Cell.Type"    "Brain.Region"          "Gene.Symbol..Human."  
[5] "Gene.Symbol..Mouse."   "Species"               "Age"                   "Statistical.Criterion"
[9] "Specificity"           "Comparison"            "Platform"              "Citation"             
[13] "Tag"                   "CellType_Primary" 

table(CellTypeSpecificGenes_Master3$CellType_Primary)

#Note: I had to reverse the order of the next 3 lines to properly remove NAs - I should change that before attempting to release the code in a way that other people will use.

colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"
colnames(CellTypeSpecificGenes_Master3)[5]<-"GeneSymbol_Mouse"

sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human))
#[1] 364
sum(is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse))
#[1] 7

#This line depends on the species being used:
CellTypeSpecificGenes_Master3NoNA<-CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse)==F,]

CellTypeSpecificGenes_Master3NoNA[,4]<-as.character(CellTypeSpecificGenes_Master3NoNA[,4])
CellTypeSpecificGenes_Master3NoNA[,5]<-as.character(CellTypeSpecificGenes_Master3NoNA[,5])


str(CellTypeSpecificGenes_Master3NoNA)

#Note to self - if I try to make a generalizable version of this code, I will need to change it so that it references the rodent gene symbols when  using rodent data and the human gene symbols when using human data.

#Note: some data will need to be averaged by gene symbol at this point.  This typically is not true of processed RNASeq data
sum(unique(row.names(ZscoreZhang))==F)
#[1] 0

#New version (after averaging by gene symbol if necessary):
temp<-data.frame(row.names(ZscoreZhang), ZscoreZhang, stringsAsFactors=F)

#This code depends on the species of the original dataset
#colnames(temp)[1]<-"GeneSymbol_Human"
colnames(temp)[1]<-"GeneSymbol_Mouse"

sum(is.na(temp[,1]))
#[1] 0

#This code also depends on the species in the original dataset:
#If Human: sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,4])
sum(temp[,1] %in% CellTypeSpecificGenes_Master3[,5])
#[1] 2513

#If human: sum(CellTypeSpecificGenes_Master3[,4]  %in%  temp[,1])
sum(CellTypeSpecificGenes_Master3[,5]  %in%  temp[,1])
# [1] 2914

#Note: NAs were causing a serious problem with this join function.  Fixed now. :)
library(plyr)
#If human: ZscoreZhang_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Human", type="inner")
ZscoreZhang_Expression_CellType<-join(CellTypeSpecificGenes_Master3, temp, by="GeneSymbol_Mouse", type="inner")
dim(ZscoreZhang_Expression_CellType)
# [1] 2914   31
#It is making all possible combinations - some of the cell type specific genes are found in more than one index.

write.csv(ZscoreZhang_Expression_CellType, "ZscoreZhang_Expression_CellType.csv")

###############################################

AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreZhang_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreZhang_Expression_CellType)-14))

row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))
colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(temp)[-1]


for(i in c(15:ncol(ZscoreZhang_Expression_CellType))){
  AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreZhang_Expression_CellType[,i], ZscoreZhang_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
}

head(AVE_Expression_CellType_Primary_bySample)

png("CorrMatrixCellTypeVsCellType_HeatMap.png")
heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
dev.off()

CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))

write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
#Huh - all correlations are positive. Perhaps because some samples simply have less reads or more artifacts?  

#let's see some examples:

png("AstrocyteByEndothelial.png")
plot(AVE_Expression_CellType_Primary_bySample[1,]~AVE_Expression_CellType_Primary_bySample[2,])
dev.off()
#I think the sample size per cell type maybe too small to properly run these analyses 

AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreZhang_Expression_CellType$Tag))), ncol=ncol(ZscoreZhang_Expression_CellType)-14)
row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreZhang_Expression_CellType$Tag))
colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(temp)[-1]

for(i in c(15:ncol(ZscoreZhang_Expression_CellType))){
  AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreZhang_Expression_CellType[,i], ZscoreZhang_Expression_CellType$Tag, mean)
}

head(AVE_Expression_CellType_Tag_bySample)

png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))

write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")

################################################

#Alright, so part of the trouble here is that there hasn't been any removal of overlapping probes yet, and we aren't averaging by tag. Let's go ahead and do that.


#Making a storage matrix to store information about overlap between primary indices:
CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(ZscoreZhang_Expression_CellType$CellType_Primary)), length(table(ZscoreZhang_Expression_CellType$CellType_Primary)) )

colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))
row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreZhang_Expression_CellType$CellType_Primary))

#Quantifying overlap between primary cell type indices:
for(i in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
  for(j in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
    
    CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4]%in%ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[j]), 4])/length(ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4])
    
  }
}

write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")



#What happens if we eliminate overlap between primary categories and then make master indices:

dim(ZscoreZhang_Expression_CellType)
# [1] 2914   31

#Making an empty first row for the storage matrix:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(ZscoreZhang_Expression_CellType[1,])))
colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)<-colnames(ZscoreZhang_Expression_CellType)

for(i in 1: length(table(ZscoreZhang_Expression_CellType$CellType_Primary))){
  
  #Choosing all data for a particular primary cell type:
  TempCurrentIndexAllInfo<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), ] 
  
  #All of the gene symbols within the current primary cell type:
  TempCurrentIndex<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary==names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[i]), 4] 
  
  #All of the gene symbols within all other primary cell types:
  TempAllOtherIndices<-ZscoreZhang_Expression_CellType[ZscoreZhang_Expression_CellType$CellType_Primary%in%names(table(ZscoreZhang_Expression_CellType$CellType_Primary)[-i]), 4]
  
  #Grabs only rows of data with gene symbols not found in other primary cell type indices:
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-rbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])
  
}

dim(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)
# [1] 2435   31

#removing that one dummy row:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap<-ZscoreZhang_Expression_CellType_NoPrimaryOverlap[-1,]

dim(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)
# [1] 2434   31

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap.csv")


CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary)

write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap.csv")

ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

#Old version of code:
# for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
# ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary, mean)
# }


#I went back and changed this so that it averaged by tag first, then by primary cell category

temp<-data.frame(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
dim(temp)
# [1] 2434    2

CellTypePrimaryVsTag<-unique(temp)
dim(CellTypePrimaryVsTag)
# [1] 38  2

colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
head(CellTypePrimaryVsTag)
#                                    Tag CellType_Primary
# 1     Astrocyte_All_Zhang_PNAS_2015        Astrocyte
# 21     Astrocyte_All_Cahoy_JNeuro_2008        Astrocyte
# 72     Astrocyte_All_Zhang_JNeuro_2014        Astrocyte
# 102      Astrocyte_All_Doyle_Cell_2008        Astrocyte
# 118  Astrocyte_All_Zeisel_Science_2015        Astrocyte
# 309 Endothelial_All_Zhang_PNAS_2015      Endothelial

ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}

head(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag.csv")


#Making histograms for each cell type tag:

for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  png(paste("Histogram_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""))
  hist(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], breaks=16, col=i)
  dev.off()
}

#Plotting each cell type tag vs. Actual cell type:
for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  png(paste("CellTypeVs_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""), height=500, width=1000)
  plot(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,]~as.factor(ExperimentalCellType), main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ylab=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], xlab="Actual Cell Type", col=i)
  dev.off()
}




#Plotting all the cell type tags for each cell type:
ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType<-matrix(0, nrow=nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag), ncol=length(unique(ExperimentalCellType)))
row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)
colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)<-names(table(ExperimentalCellType))

for(i in 1: nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[i,]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], ExperimentalCellType, mean)
}

write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType.csv")

ColorsForCellTypes<-as.character(CellTypePrimaryVsTag[,2])
ColorsForCellTypes[ColorsForCellTypes=="Astrocyte"]<-"lavender"
ColorsForCellTypes[ColorsForCellTypes=="Endothelial"]<-"orange"
ColorsForCellTypes[ColorsForCellTypes=="Microglia"]<-"darkolivegreen3"
ColorsForCellTypes[ColorsForCellTypes=="Mural"]<-"yellow"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_All"]<-"darkviolet"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Interneuron"]<-"blue"
ColorsForCellTypes[ColorsForCellTypes=="Neuron_Projection"]<-"red"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte"]<-"lightpink"
ColorsForCellTypes[ColorsForCellTypes=="Oligodendrocyte_Immature"]<-"ivory2"
ColorsForCellTypes[ColorsForCellTypes=="RBC"]<-"indianred4"

#There is a mysterious legend that keeps popping up... that is just the *original* margin parameters. Um. Yeah. Weird, right?
for(i in 1:length(unique(ExperimentalCellType))){
  png(paste("CellTypeTagVs_", colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ".png", sep=""), height=800, width=1500)
  barplot(height=ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType[,i], width=1, space = NULL, main=colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType)[i], ylab="Average Cell Type Index (Mean Z-score for All Genes)", xlab=NULL, names.arg=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType), las=2, par(mar=c(25,6,4,2)+0.1), cex.lab=1.5, col=ColorsForCellTypes)
  dev.off()
}

temp<-apply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,c(9:10,16:17)], 1, mean)

png(paste("CellTypeTagVs_EnrichedForWhiteMatter2", ".png", sep=""), height=800, width=1500)
barplot(height=temp, width=1, space = NULL, main="Whole Brain Enriched for WM (WB=2, WM=2)", ylab="Average Cell Type Index (Mean Z-score for All Genes)", xlab=NULL, names.arg=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType), las=2, par(mar=c(25,6,4,2)+0.1), cex.lab=1.5, col=ColorsForCellTypes)
dev.off()

temp<-apply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,c(9,16:17)], 1, mean)

png(paste("CellTypeTagVs_EnrichedForWhiteMatter1", ".png", sep=""), height=800, width=1500)
barplot(height=temp, width=1, space = NULL, main="Whole Brain Enriched for WM (WB=2, WM=1)", ylab="Average Cell Type Index (Mean Z-score for All Genes)", xlab=NULL, names.arg=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag_byActualCellType), las=2, par(mar=c(25,6,4,2)+0.1), cex.lab=1.5, col=ColorsForCellTypes)
dev.off()

#determining how tightly we can track in silico ratios:

AstrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(AstrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
AstrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(17:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  AstrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(AstrocyteMixtures)

plot(AstrocyteMixtures[1,]~AstrocyteMixtureRatios, pch=18)
points(AstrocyteMixtures[2,]~AstrocyteMixtureRatios, col=2, pch=18)
points(AstrocyteMixtures[3,]~AstrocyteMixtureRatios, col=3, pch=18)
points(AstrocyteMixtures[4,]~AstrocyteMixtureRatios, col=4, pch=18)
points(AstrocyteMixtures[5,]~AstrocyteMixtureRatios, col=5, pch=18)

png("AstrocyteMixtureRatiosVsIndex_NoZhang.png")
plot(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="AstrocyteIndex")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()


#Version 2: Using whole brain as the background instead of a random sampling of cells
AstrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(AstrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
AstrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(29:31), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  AstrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(AstrocyteMixtures)

plot(AstrocyteMixtures[1,]~AstrocyteMixtureRatios, pch=18)
points(AstrocyteMixtures[2,]~AstrocyteMixtureRatios, col=2, pch=18)
points(AstrocyteMixtures[3,]~AstrocyteMixtureRatios, col=3, pch=18)
points(AstrocyteMixtures[4,]~AstrocyteMixtureRatios, col=4, pch=18)
points(AstrocyteMixtures[5,]~AstrocyteMixtureRatios, col=5, pch=18)

png("AstrocyteMixtureRatiosVsIndex_NoZhang_WB.png")
plot(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios, pch=18, ylab="AstrocyteIndex")
BestFitLine<-lm(apply(AstrocyteMixtures[c(1:4),], 2, mean)~AstrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()



OligodendrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(OligodendrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
OligodendrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(23:24), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:22, 25:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  OligodendrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(OligodendrocyteMixtures)

png("OligodendrocyteMixtureRatiosVsIndex_NoZhang.png")
plot(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="OligodendrocyteIndex")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()


OligodendrocyteMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(OligodendrocyteMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
OligodendrocyteMixtureRatios<-c(rep(15, 5), rep(20, 5), rep(25, 5), rep(30,5), rep(35,5), rep(40,5))

for(i in c(15, 20, 25, 30, 35, 40)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(23:24), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(29:31), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  OligodendrocyteMixtures[,c(i:(i+4))-14]<-temp3
}

head(OligodendrocyteMixtures)

png("OligodendrocyteMixtureRatiosVsIndex_NoZhang_WB.png")
plot(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios, pch=18, ylab="OligodendrocyteIndex")
BestFitLine<-lm(apply(OligodendrocyteMixtures[c(29:33),], 2, mean)~OligodendrocyteMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*30), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*27.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*25), b=0, col="red")
dev.off()

#I'm going to stop doing the whole brain mixtures - i think there isn't a reasonable amount of variability there, since there were only 3 WB samples.  I general, that is a definite problem with this entire analysis.


MicrogliaMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(MicrogliaMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
MicrogliaMixtureRatios<-c(rep(5, 5), rep(10, 5), rep(15, 5), rep(20,5), rep(25,5), rep(30,5))

for(i in c(5, 10, 15, 20, 25, 30)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(25:26), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:24, 27:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  MicrogliaMixtures[,c(i:(i+4))-4]<-temp3
}

head(MicrogliaMixtures)

plot(MicrogliaMixtures[1,]~MicrogliaMixtureRatios, pch=18)
points(MicrogliaMixtures[2,]~MicrogliaMixtureRatios, col=2, pch=18)
points(MicrogliaMixtures[3,]~MicrogliaMixtureRatios, col=3, pch=18)
points(MicrogliaMixtures[4,]~MicrogliaMixtureRatios, col=4, pch=18)
points(MicrogliaMixtures[5,]~MicrogliaMixtureRatios, col=5, pch=18)

png("MicrogliaMixtureRatiosVsIndex_NoZhang.png")
plot(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios, pch=18, ylab="MicrogliaIndex")
BestFitLine<-lm(apply(MicrogliaMixtures[c(10:11),], 2, mean)~MicrogliaMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*15), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*17.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*20), b=0, col="red")
dev.off()

EndothelialMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(EndothelialMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
EndothelialMixtureRatios<-c(rep(1, 5), rep(6, 5), rep(11, 5), rep(16,5), rep(21,5), rep(26,5))

for(i in c(1, 6, 11, 16, 21, 26)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(27:28), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:26), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  EndothelialMixtures[,c(i:(i+4))]<-temp3
}

head(EndothelialMixtures)

plot(EndothelialMixtures[1,]~EndothelialMixtureRatios, pch=18)
points(EndothelialMixtures[2,]~EndothelialMixtureRatios, col=2, pch=18)
points(EndothelialMixtures[3,]~EndothelialMixtureRatios, col=3, pch=18)
points(EndothelialMixtures[4,]~EndothelialMixtureRatios, col=4, pch=18)
points(EndothelialMixtures[5,]~EndothelialMixtureRatios, col=5, pch=18)

png("EndothelialMixtureRatiosVsIndex_NoZhang.png")
plot(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios, pch=18, ylab="EndothelialIndex")
BestFitLine<-lm(apply(EndothelialMixtures[c(6:8),], 2, mean)~EndothelialMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*5), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*7.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*10), b=0, col="red")
dev.off()

NeuronMixtures<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 30)
row.names(NeuronMixtures)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))
NeuronMixtureRatios<-c(rep(25, 5), rep(30, 5), rep(35, 5), rep(40,5), rep(45,5), rep(50,5))

for(i in c(25, 30, 35, 40, 45, 50)){
  
  temp3<-matrix(0, length(names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))), 5)  
  
  for(j in 1:5){
    temp<-cbind(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(17:18), i, replace=T)],
                ZscoreZhang_Expression_CellType_NoPrimaryOverlap[, sample(c(15:16, 19:28), (100-i), replace=T)])
    temp2<-apply(temp, 1, mean)
    temp3[,j]<-tapply(temp2, ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
  }
  
  NeuronMixtures[,c(i:(i+4))-24]<-temp3
}

head(NeuronMixtures)

plot(NeuronMixtures[1,]~NeuronMixtureRatios, pch=18)
points(NeuronMixtures[2,]~NeuronMixtureRatios, col=2, pch=18)
points(NeuronMixtures[3,]~NeuronMixtureRatios, col=3, pch=18)
points(NeuronMixtures[4,]~NeuronMixtureRatios, col=4, pch=18)
points(NeuronMixtures[5,]~NeuronMixtureRatios, col=5, pch=18)

png("NeuronMixtureRatiosVsIndex_NoZhang.png")
plot(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios, pch=18, ylab="NeuronIndex")
BestFitLine<-lm(apply(NeuronMixtures[c(16:17),], 2, mean)~NeuronMixtureRatios)
abline(BestFitLine)
abline(a=(BestFitLine$coefficients[1]+summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")
abline(a=(BestFitLine$coefficients[1]-summary.lm(BestFitLine)[[6]]), b=BestFitLine$coefficients[2], col="grey")

abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*35), b=0, col="red")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*37.5), b=0, col="green")
abline(a=BestFitLine$coefficients[1]+(BestFitLine$coefficients[2]*40), b=0, col="red")
dev.off()


ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,])-14))

row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag))

colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]

for(i in c(15:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,i], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
}


################################

png("Heatmap_CellType_NoPrimaryOverlap_MeanTag.png", height=1000, width=1000)
heatmap(cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellType_NoPrimaryOverlap_MeanTag_CorrMatrix<-cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F]))

head(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix)

write.csv(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix, "CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv")
###########

temp2<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap[,15], ZscoreZhang_Expression_CellType_NoPrimaryOverlap$Tag, mean)
Tag<-names(temp2)

CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")


for(i in c(1:length(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
  ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
}

head(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)


write.csv(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean, "ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean.csv")

#Making histograms for each primary cell type:

for(i in 1:nrow(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)){
  png(paste("Histogram_", row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)[i], ".png", sep=""))
  hist(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[i,], main=row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)[i], breaks=16, col=i)
  dev.off()
}

#Interesting - it is easy to say what is *not* the cell type of interest, but the values for what could be the cell type of interest range greatly. I'm guessing that this is partially a property of the skewed variability and signal values in the data itself, but I'm not sure. I wonder if it correlates at all with the read qc stats for the samples.


is.numeric(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)

png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png", height=1000, width=1000)
heatmap(cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])), cex.lab=0.3, margins = c(20, 20))
dev.off()

CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))

head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)

write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")

##################################################

#O.k., note that this analysis so far is somewhat self-referential, since genes that Zhang found to be cell type specific are included. 

#Also, without having an external reference, it is hard to interpret them. Many of the distributions sort of fit expectation... but that's not really enough to make an argument.

#So that leaves us with either Cahoy or Zhang. If we use either of them, we don't really have a good way to apply standard deconvolution techniques. hmmm... And Cahoy doesn't include the minor cell types.
