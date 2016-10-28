cellTypeFunction <- function(userInput, dataColumns, geneColumn, species){
  
  #required packages
  #install.packages("gdata")
  #library(gdata)
  #install.packages("fields")
  #library(fields)
  
  #library(stats)
  #install.packages("car")
  #library(car)
  
  #FIXME
  #not availible for 3.3.1 
  #install.packages("affy")
  #library(affy)
  #install.packages("preprocessCore")
  #library(preprocessCore)
  #install.packages("multtest")
  #library(multtest)
  
  #DATASET READ IN
    #userInput <- read.table(file = fileName, header = T, sep = ",", stringsAsFactors = F)
    #userInput <- read.table(file = "JoinedZhang_AsNumMatrix_Log2.csv", header = T, sep = ",", stringsAsFactors = F)
    GeneNamesForJoinedInput <- userInput[,1]
    GeneNamesForJoinedInput <- as.matrix(GeneNamesForJoinedInput)
    
  #########
    #chagning values based on user Input
    #if (species == "Mouse" || species == "mouse"){
    #  CellTypeGeneColumn = 5
    #}
    #else{
    #  CellTypeGeneColumn = 4
    #}
    
    #function begins
    
    ##################################################
    #correlation matrices
    TempJoinedInput_AsNum <- userInput[,2:18]
    #temp<-cor(TempJoinedInput_AsNum)
    #row.names(temp)<-colnames(userInput)
    
    #png("09 Sample Sample Correlations Heatmap.png")
    #  heatmap(temp, main="Visualizing correlations between entire samples", xlab="Red=Less correlated, Light yellow=Highly correlated")
    #dev.off()
    #Note that the heatmap can be tailored to focus on a certain level of correlation by using the command zlim=c(lower limit, upper limit)
    
    #Visualize the sample-sample correlations using a boxplot:
    #png("09 Boxplot Sample Sample Correlations.png", width=2000, height=300)
    #  boxplot(data.frame(cor(TempJoinedInput_AsNum)), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of sample-sample correlations", xlab="Subject", ylab="Sample-Sample Correlations")
    #  Median10thQuantile<-median(apply((cor(TempJoinedInput_AsNum)), 1, quantile, 0.1))
    #  MedianQuantile<-median(apply((cor(TempJoinedInput_AsNum)), 1, quantile, 0.5))
    #  abline(a=Median10thQuantile, b=0, col=2)
    #  abline(a=MedianQuantile, b=0, col=3)
    #  mtext(paste("Median Sample-Sample Correlation=", round(MedianQuantile, digits=3), sep=" ")) 
    #dev.off()
    
    
    
    ##################################################
    #Test stuff
    
    temp3<-apply(TempJoinedInput_AsNum, 1, max)
        sum(temp3<2)
        #9787
    
        sum(temp3<3)
        #11889
    
        sum(temp3<3)/length(temp3)
        #.529 or 53%
    
    ##################################################
    JoinedInput_StDev<-apply(TempJoinedInput_AsNum, 1, sd) 
        sum(JoinedInput_StDev==0)
        #5314
    
    
    JoinedInput_AsNumMatrix_Log2_NoSD0<-TempJoinedInput_AsNum[JoinedInput_StDev>0,]
        temp<-GeneNamesForJoinedInput
        GeneNamesForJoinedInput_NoSD0<-temp[JoinedInput_StDev>0]
        GeneNamesForJoinedInput_NoSD0 <- as.matrix(GeneNamesForJoinedInput_NoSD0)
    ZscoreInput<-t(scale(t(JoinedInput_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
        write.csv(ZscoreInput, "ZscoreInput.csv")

        sum(is.na(ZscoreInput))
        # ZERO WOOOOOOOOOO
    
    #############################################################
    
    CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)
    colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"
    colnames(CellTypeSpecificGenes_Master3)[5]<-"GeneSymbol_Mouse"
    #removing na values
    
    
    
    #############################################################
    #joining celltype to zscore data
    tempForJoin <- data.frame(GeneNamesForJoinedInput_NoSD0, ZscoreInput, stringsAsFactors=F)
    library(plyr)
    
    #########
    #chagning values based on user Input
    if (species == "Mouse" || species == "mouse"){
      CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse)==F,]
      colnames(CellTypeSpecificGenes_Master3NoNA)[5] <- "GeneNamesForJoinedInput_NoSD0"
      ZscoreInput_Expression_CellType<<-join(CellTypeSpecificGenes_Master3NoNA, tempForJoin, by="GeneNamesForJoinedInput_NoSD0", type="inner")
      write.csv(ZscoreInput_Expression_CellType, "ZscoreInput_Expression_CellType.csv")
    }
    else if(species == "Human" || species == "human"){
      CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Human)==F,]
      colnames(CellTypeSpecificGenes_Master3NoNA)[4] <- "GeneNamesForJoinedInput_NoSD0"
      ZscoreInput_Expression_CellType<-join(CellTypeSpecificGenes_Master3NoNA, tempForJoin, by="GeneNamesForJoinedInput_NoSD0", type="inner")
      write.csv(ZscoreInput_Expression_CellType, "ZscoreInput_Expression_CellType.csv")
    }

    ####################################
    #///////////////////////////////////
    ####################################
    #do these need to get output or are they testing?
    #CELLTYPE PRIMARY BY SAMPLE
    AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreInput_Expression_CellType)-14))
    row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
    colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(tempForJoin)[-1]
    
    for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
      AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
    }
    
    png("CorrMatrixCellTypeVsCellType_HeatMap.png")
    heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
    dev.off()
    
    CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))
    
    write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
    #Huh - all correlations are positive. Perhaps because some samples simply have less reads or more artifacts? 

    #####################################
    # CELL TYPE TAG BY SAMPLE
    AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$Tag))), ncol=ncol(ZscoreInput_Expression_CellType)-14)
    row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreInput_Expression_CellType$Tag))
    colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(tempForJoin)[-1]
    
    for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
      AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$Tag, mean)
    }
    
    png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
    heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
    dev.off()
    
    CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))
    
    write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")
    
    #############################################
    #////////////////////////////////////////////
    #############################################
    
    #Alright, so part of the trouble here is that there hasn't been any removal of overlapping probes yet, and we aren't averaging by tag. Let's go ahead and do that.
    
    
    #Making a storage matrix to store information about overlap between primary indices:
    CellTypeSpecificGenes_Master3_Overlap<-matrix(0, length(table(ZscoreInput_Expression_CellType$CellType_Primary)), length(table(ZscoreInput_Expression_CellType$CellType_Primary)) )
    
    colnames(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
    row.names(CellTypeSpecificGenes_Master3_Overlap)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
    
    #Quantifying overlap between primary cell type indices:
    for(i in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
        for(j in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
            
            CellTypeSpecificGenes_Master3_Overlap[i,j]<-sum(ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4]%in%ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[j]), 4])/length(ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4])
            
        }
    }
    
    write.csv(CellTypeSpecificGenes_Master3_Overlap, "CellTypeSpecificGenes_Master3_Overlap.csv")
    
    
    
    #############################################
    #What happens if we eliminate overlap between primary categories and then make master indices:
    
    dim(Zscorenput_Expression_CellType)
    # [1] 2914   31
    
    #Making an empty first row for the storage matrix:
    #
    ZscoreInput_Expression_CellType_NoPrimaryOverlap<-matrix(0, 1, (length(ZscoreInput_Expression_CellType[1,])))
    colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)<-colnames(ZscoreInput_Expression_CellType)
    
    for(i in 1: length(table(ZscoreInput_Expression_CellType$CellType_Primary))){
        
        #Choosing all data for a particular primary cell type:
        TempCurrentIndexAllInfo<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), ] 
        
        #All of the gene symbols within the current primary cell type:
        TempCurrentIndex<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary==names(table(ZscoreInput_Expression_CellType$CellType_Primary)[i]), 4] 
        
        #All of the gene symbols within all other primary cell types:
        TempAllOtherIndices<-ZscoreInput_Expression_CellType[ZscoreInput_Expression_CellType$CellType_Primary%in%names(table(ZscoreInput_Expression_CellType$CellType_Primary)[-i]), 4]
        
        #Grabs only rows of data with gene symbols not found in other primary cell type indices:
        ZscoreInput_Expression_CellType_NoPrimaryOverlap<-rbind(ZscoreInput_Expression_CellType_NoPrimaryOverlap, TempCurrentIndexAllInfo[(TempCurrentIndex%in%TempAllOtherIndices)==F,])
        
    }
    
    dim(ZscoreInput_Expression_CellType_NoPrimaryOverlap)
    # [1] 2435   31
    
    #removing that one dummy row:
    ZscoreInput_Expression_CellType_NoPrimaryOverlap<-ZscoreInput_Expression_CellType_NoPrimaryOverlap[-1,]
    
    dim(ZscoreInput_Expression_CellType_NoPrimaryOverlap)
    # [1] 2434   31
    
    write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap, "ZscoreInput_Expression_CellType_NoPrimaryOverlap.csv")
    
    ######################################################################
    #HALP, WHAT DO?
    #CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap<-table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary)
    
    
    #write.csv(CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap, "CellTypeSpecificGenes_Master3_PrimaryTable_NoPrimaryOverlap.csv")
    
    #ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,])-14))
    
    #row.names(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary))
    
    #colnames(ZscoreZhang_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]
    
    #end help what do
    # Megan has an "old version" forloop commented out in the Original doc.
    #is this function not used?
    ######################################
    
    temp<-data.frame(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary) 
    dim(temp)
    # [1] 2434    2
    
    CellTypePrimaryVsTag<-unique(temp)
    dim(CellTypePrimaryVsTag)
    # [1] 38  2
    
    colnames(CellTypePrimaryVsTag)<-c("Tag","CellType_Primary")
    head(CellTypePrimaryVsTag)
    
    #####################################################
    ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag<-matrix(0, nrow=length(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag)), ncol=(length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,])-14))
    
    row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)<-names(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag))
    
    colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)<-colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]
    
    for(i in c(15:length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,]))){
        ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,i-14]<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap[,i], ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, mean)
    }
    
    head(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)
    
    write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag, "Zscore_Expression_CellType_NoPrimaryOverlap_MeanTag.csv")
    
    #Making histograms for each cell type tag:
    dir.create("Cell Type Histograms")
    setwd("Cell Type Histograms")
    for(i in 1:nrow(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)){
        png(paste("Histogram_", row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], ".png", sep=""))
        hist(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[i,], main=row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag)[i], breaks=16, col=i)
        dev.off()
    }
    setwd("..")
    getwd()
    
    ############################33
    
    png("Heatmap_CellType_NoPrimaryOverlap_MeanTag.png", height=1000, width=1000)
    heatmap(cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F])), cex.lab=0.3, margins = c(20, 20))
    dev.off()
    
    CellType_NoPrimaryOverlap_MeanTag_CorrMatrix<-cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,])==F]))
    
    head(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix)
    
    write.csv(CellType_NoPrimaryOverlap_MeanTag_CorrMatrix, "CellType_NoPrimaryOverlap_MeanTag_CorrMatrix.csv")
    
    
    
    ###############################################
    
    temp2<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap[,15], ZscoreInput_Expression_CellType_NoPrimaryOverlap$Tag, mean)
    Tag<-names(temp2)
    
    CellTypePrimaryVsTag2<-join(as.data.frame(Tag), as.data.frame(CellTypePrimaryVsTag), by="Tag")
    
    ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean<-matrix(0, nrow=length(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary)), ncol=(length(ZscoreInput_Expression_CellType_NoPrimaryOverlap[1,])-14))
    
    row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)<-names(table(ZscoreInput_Expression_CellType_NoPrimaryOverlap$CellType_Primary))
    
    colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)<-colnames(ZscoreInput_Expression_CellType_NoPrimaryOverlap)[-c(1:14)]
    
    for(i in c(1:length(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[1,]))){
        ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,i]<-tapply(ZscoreInput_Expression_CellType_NoPrimaryOverlap_MeanTag[,i], CellTypePrimaryVsTag2[,2], mean)
    }
    
    head(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)
    
    write.csv(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean, "Zscore_Expression_CellType_NoPrimaryOverlap_Mean.csv")
    
    
    ##############################################
    #Making histograms for each primary cell type:
    dir.create("Primary cell type Histogram")
    setwd("Primary cell type Histogram")
    for(i in 1:nrow(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)){
        png(paste("Histogram_", row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)[i], ".png", sep=""))
        hist(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[i,], main=row.names(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)[i], breaks=16, col=i)
        dev.off()
    }
    
    #Interesting - it is easy to say what is *not* the cell type of interest, but the values for what could be the cell type of interest range greatly. I'm guessing that this is partially a property of the skewed variability and signal values in the data itself, but I'm not sure. I wonder if it correlates at all with the read qc stats for the samples.
    
    
    is.numeric(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean)
    
    png("Heatmap_CorMatrixPrimaryCellsNoOverlap.png", height=1000, width=1000)
    heatmap(cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F])), cex.lab=0.3, margins = c(20, 20))
    dev.off()
    
    CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix<-cor(t(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[,is.na(ZscoreInput_Expression_CellType_NoPrimaryOverlap_Mean[1,])==F]))
    
    head(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix)
    
    write.csv(CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix, "CellTypeIndexNoNA3_NormBest_NoPrimaryOverlap_CorrMatrix.csv")
    
    ##################################################
    
  return(CellTypeSpecificGenes_Master3NoNA)
}

Zhangdoc <- read.table(file = "JoinedZhang_AsNumMatrix_Log2.csv", header = T, sep = ",", stringsAsFactors = F)

dump("cellTypeFunction", file = "CellTypeFunction.R")
source("cellTypeFunction.R")

t<-cellTypeFunction(userInput = Zhangdoc, dataColumns = 2:18, geneColumn = 1, species = "Human")

t <- as.data.frame(t)

levels(CellTypeSpecificGenes_Master3$CellType_Primary)
