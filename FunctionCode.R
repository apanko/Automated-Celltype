cellTypeFunction <- function(fileName, dataColumns, geneColumn, Species){
  
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
    userInput <- read.table(file = fileName, header = T, sep = ",", stringsAsFactors = F)
    #userInput <- read.table(file = "JoinedZhang_AsNumMatrix_Log2.csv", header = T, sep = ",", stringsAsFactors = F)
    GeneNamesForJoinedInput <- userInput[,geneColumn]
    GeneNamesForJoinedInput <- as.matrix(GeneNamesForJoinedInput)
    
  #########
    #chagning values based on user Input
    if (species == "Mouse" || species == "mouse"){
      CellTypeGeneColumn = 5
    }
    else{
      CellTypeGeneColumn = 4
    }
    
    #function begins
    
    ##################################################
    #correlation matrices
    #TempJoinedInput_AsNum <- userInput[,2:18]
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
    
        sum(temp3<3)/length(temp)
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
        ZscoreInput <- as.table(ZscoreInput)
        sum(is.na(ZscoreInput))
        # ZERO WOOOOOOOOOO
    
    #############################################################
    
    CellTypeSpecificGenes_Master3<-read.csv("CellTypeSpecificGenes_Master3.csv", header=T)
    colnames(CellTypeSpecificGenes_Master3)[4]<-"GeneSymbol_Human"
    colnames(CellTypeSpecificGenes_Master3)[5]<-"GeneSymbol_Mouse"
    #removing na values
    CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3[is.na(CellTypeSpecificGenes_Master3$GeneSymbol_Mouse)==F,]
    CellTypeSpecificGenes_Master3NoNA <- CellTypeSpecificGenes_Master3NoNA[is.na(CellTypeSpecificGenes_Master3NoNA$GeneSymbol_Human)==F,]
    
    #############################################################
    #joining celltype to zscore data
    temp<-data.frame(GeneNamesForJoinedInput_NoSD0, ZscoreInput, stringsAsFactors=F)
    
    
    #########
    #chagning values based on user Input
    if (species == "Mouse" || species == "mouse"){
      colnames(CellTypeSpecificGenes_Master3NoNA)[5] <- "GeneNamesForJoinedInput_NoSD0"
    }
    else{
      colnames(CellTypeSpecificGenes_Master3NoNA)[4] <- "GeneNamesForJoinedInput_NoSD0"
      
    }

    
  return(CellTypeSpecificGenes_Master3NoNA)
}

dump("cellTypeFunction", file = "CellTypeFunction.R")
source("cellTypeFunction.R")

t<-cellTypeFunction(fileName = "JoinedZhang_AsNumMatrix_Log2.csv", dataColumns = 2:18, geneColumn = 1, Species = "Mouse")

t <- as.data.frame(t)

