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

    ################################3
    
    AVE_Expression_CellType_Primary_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$CellType_Primary))), ncol=(ncol(ZscoreInput_Expression_CellType)-14))
    row.names(AVE_Expression_CellType_Primary_bySample)<-names(table(ZscoreInput_Expression_CellType$CellType_Primary))
    colnames(AVE_Expression_CellType_Primary_bySample)<-colnames(temp)[-1]
    
    for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
      AVE_Expression_CellType_Primary_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$CellType_Primary, function(y) mean(y, na.rm=T))
    }
    
    png("CorrMatrixCellTypeVsCellType_HeatMap.png")
    heatmap(cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F])), margins = c(15, 15), cex.lab=0.5)
    dev.off()
    
    CorrelationMatrixCellTypeVsCellType<-cor(t(AVE_Expression_CellType_Primary_bySample[,is.na(AVE_Expression_CellType_Primary_bySample[1,])==F]))
    
    write.csv(CorrelationMatrixCellTypeVsCellType, "CorrelationMatrixCellTypeVsCellType.csv")
    #Huh - all correlations are positive. Perhaps because some samples simply have less reads or more artifacts? 
    
    AVE_Expression_CellType_Tag_bySample<-matrix(NA, nrow=length(names(table(ZscoreInput_Expression_CellType$Tag))), ncol=ncol(ZscoreInput_Expression_CellType)-14)
    row.names(AVE_Expression_CellType_Tag_bySample)<-names(table(ZscoreInput_Expression_CellType$Tag))
    colnames(AVE_Expression_CellType_Tag_bySample)<-colnames(temp)[-1]
    
    for(i in c(15:ncol(ZscoreInput_Expression_CellType))){
      AVE_Expression_CellType_Tag_bySample[,(i-14)]<-tapply(ZscoreInput_Expression_CellType[,i], ZscoreInput_Expression_CellType$Tag, mean)
    }
    
    png("CorrMatrixCellIndexVsCellIndex_HeatMap.png", width=1000, height=1000)
    heatmap(cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F])), cex.lab=0.3, margins = c(20, 20))
    dev.off()
    
    CorrelationMatrixCellIndexVsCellIndex<-cor(t(AVE_Expression_CellType_Tag_bySample[,is.na(AVE_Expression_CellType_Tag_bySample[1,])==F]))
    
    write.csv(CorrelationMatrixCellIndexVsCellIndex, "CorrelationMatrixCellIndexVsCellIndex.csv")
    
  return(CellTypeSpecificGenes_Master3NoNA)
}

Zhangdoc <- read.table(file = "JoinedZhang_AsNumMatrix_Log2.csv", header = T, sep = ",", stringsAsFactors = F)

dump("cellTypeFunction", file = "CellTypeFunction.R")
source("cellTypeFunction.R")

t<-cellTypeFunction(userInput = Zhangdoc, dataColumns = 2:18, geneColumn = 1, species = "Human")

t <- as.data.frame(t)


