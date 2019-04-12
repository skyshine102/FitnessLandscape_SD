######################################################
##                                                  ##
##  Author: Antony Kuo, National Taiwan University  ##
##  Last update: 2019/04/12                         ##
##                                                  ##
######################################################


### Substitution Dataframe ====
## (i,j): 0A, 1C, 2G, 3U
SubDF <- function(Exp,site,i,j){
  if(i == j){return(NA)}
  n <- length(Exp)
  index1 <- which(0:(n-1)%/%4^(9-site)%%4 == i)
  index2 <- which(0:(n-1)%/%4^(9-site)%%4 == j)
  DF <- data.frame(Exp[index1],Exp[index2] - Exp[index1])
  colnames(DF) <- c("x","y")
  return(DF)
}

### Boxplot Dataframe ====
BoxDF <- function(Exp,DF){
  ## Clear NA
  DF <- DF[which(!is.na(DF$x)),]
  DF <- DF[which(!is.na(DF$y)),]
  
  ## Boxplot Dataframe
  n <- 20
  BDF <- data.frame(matrix(NA,n,5))
  colnames(BDF) <- c("XMean","YMean","YQ1","YMedian","YQ3")
  Max <- 3.2
  Min <- 0
  R <- Max - Min
  r <- R/n
  for(k in 0:(n-1)){
    indexBG <- which(DF$x>=k*r+Min & DF$x<=(k+1)*r+Min)
    BDF$XMean[k+1] <- mean(DF$x[indexBG])
    BDF$YMean[k+1] <- mean(DF$y[indexBG])
    BDF$YQ1[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.25, type=2))
    BDF$YMedian[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.5, type=2))
    BDF$YQ3[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.75, type=2))
  }
  return(BDF)
}

Vienna_BoxDF <- function(Exp,DF){
  ## Clear NA
  DF <- DF[which(!is.na(DF$x)),]
  DF <- DF[which(!is.na(DF$y)),]
  
  ## Boxplot Dataframe
  n <- 20
  BDF <- data.frame(matrix(NA,n,5))
  colnames(BDF) <- c("XMean","YMean","YQ1","YMedian","YQ3")
  Max <- 14
  Min <- 0
  R <- Max - Min
  r <- R/n
  for(k in 0:(n-1)){
    indexBG <- which(DF$x>=k*r+Min & DF$x<=(k+1)*r+Min)
    BDF$XMean[k+1] <- mean(DF$x[indexBG])
    BDF$YMean[k+1] <- mean(DF$y[indexBG])
    BDF$YQ1[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.25, type=2))
    BDF$YMedian[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.5, type=2))
    BDF$YQ3[k+1] <- as.numeric(quantile(DF$y[indexBG], 0.75, type=2))
  }
  
  return(BDF)
}

################################################ SubstitutionGtoC ################################################
SubstitutionGtoC <- function(f){
  Data <- read.table(file=paste0("../../data/Read Count Table/UnionTable/",f,".csv"), header=T, sep=",")
  Exp <- Data$LogMean
  
  i <- 2
  j <- 1
  
  PlotMat <- data.frame(matrix(data = NA, nrow = 20, ncol = 10))
  colnames(PlotMat) <- c("Initial Expression",paste0("Expression Change (Position",1:9,")"))
  
  BinMean <- seq(0,3.2,3.2/20)
  BinMean <- rowMeans(cbind(BinMean[1:20],BinMean[2:21]))
  PlotMat[,1] <- BinMean
  
  for(site in 1:9){
    SDF <- SubDF(Exp,site,i,j)
    BDF <- BoxDF(Exp, DF=SDF)
    
    PlotMat[,site+1] <- BDF$YMean
  }
  write.table(x = PlotMat, file = paste0("../../Result/R_Substitution/",f,"_GtoC_Mean.csv"), 
              sep = ",", quote = F, row.names = F)
}

SubstitutionGtoC("dmsC_SDR_union_count25")
SubstitutionGtoC("arti_SDR_union_count25")
SubstitutionGtoC("fepB_SDR_union_count25")

################################################ SubstitutionTotal ################################################
SubstitutionTotal <- function(f,spacing){
  Data <- read.table(file=paste0("../../data/Read Count Table/UnionTable/",f,".csv"), header=T, sep=",")
  Exp <- Data$LogMean
  
  char <- c("A","C","G","U")
  MA <- max(Exp, na.rm=T)
  MI <- min(Exp, na.rm=T)
  R <- MA - MI
  
  BaseIndex <- matrix(data = c(0,0,0,3,3,3,2,2,2,1,1,1,3,2,1,0,2,1,0,3,1,0,3,2), 
                      nrow = 12, ncol = 2)
  tiff(filename = paste0("../../Result/R_Substitution/",f,"_total.tiff"), width = 1800, height = 4800, 
       units = "px")
  par(mai = c(.1, .1, .1, .1), mfrow = c(12,9))
  
  for(k in 1:12){
    i <- BaseIndex[k,1]
    j <- BaseIndex[k,2]
    
    for(site in 1:9){
      SDF <- SubDF(Exp,site,i,j)
      BDF <- BoxDF(Exp, DF=SDF)
      
      plot(SDF$x, SDF$y, type="p", pch=".", col="gray64", xlim=c(0,3.2), ylim=c(-3.2,3.2), 
           xlab="", ylab="", main="", axes=FALSE, xaxs="i", yaxs="i")
      axis(side=1, at=seq(0, 3.2, 3.2), cex.axis=1.5, font=2, labels = F)
      axis(side=2, at=seq(-3.2, 3.2, 3.2), cex.axis=1.5, font=2, labels = F)
      box()
      abline(h=0)
      arrows(x0=BDF$XMean, y0=BDF$YQ1, x1=BDF$XMean, y1=BDF$YQ3, lwd=6, code=0, 
             col="green") # Q1-Q3
      points(BDF$XMean, BDF$YMedian, pch=16, col="blue", cex=2.5) # Median
      points(BDF$XMean, BDF$YMean, pch=16, col="red", cex=2.5) # Mean
    }
  }
  dev.off()
}

SubstitutionTotal("dmsC_SDR_union_count25",4)
SubstitutionTotal("arti_SDR_union_count25",5)
SubstitutionTotal("fepB_SDR_union_count25",6)

################################################ Vienna_SubstitutionGtoC_SeqExp ################################################
Vienna_SubstitutionTotal_Duplex <- function(f){
  Data <- read.table(file=paste0("../../data/Vienna Package Prediction/",f,"_SeqExp.csv"), header=T, sep=",")
  Exp <- -Data$LogMean
  
  i <- 2
  j <- 1
  
  PlotMat <- data.frame(matrix(data = NA, nrow = 20, ncol = 10))
  colnames(PlotMat) <- c("Initial Expression",paste0("Expression Change (Position",1:9,")"))
  
  BinMean <- seq(0,14,14/20)
  BinMean <- rowMeans(cbind(BinMean[1:20],BinMean[2:21]))
  PlotMat[,1] <- BinMean
  
  for(site in 1:9){
    SDF <- SubDF(Exp,site,i,j)
    BDF <- Vienna_BoxDF(Exp, DF=SDF)
    
    PlotMat[,site+1] <- BDF$YMean
  }
  write.table(x = PlotMat, file = paste0("../../Result/R_Substitution/",f,"_SeqExp_GtoC_Mean.csv"), 
              sep = ",", quote = F, row.names = F)
}

Vienna_SubstitutionTotal_Duplex("Vienna")

################################################ Vienna_SubstitutionTotal_SeqExp ################################################
Vienna_SubstitutionTotal_SeqExp <- function(f){
  Data <- read.table(file=paste0("../../data/Vienna Package Prediction/",f,"_SeqExp.csv"), header=T, sep=",")
  Exp <- -Data$LogMean
  
  char <- c("A","C","G","U")
  BaseIndex <- matrix(data=c(0,0,0,3,3,3,2,2,2,1,1,1,3,2,1,0,2,1,0,3,1,0,3,2), nrow=12, ncol=2)
  tiff(filename = paste0("../../Result/R_Substitution/",f,"_SeqExp_total.tiff"), width = 1800, height = 4800, units = "px")
  par(mai = c(.1, .1, .1, .1), mfrow = c(12,9))
  
  for(k in 1:12){
    i <- BaseIndex[k,1]
    j <- BaseIndex[k,2]
    
    for(site in 1:9){
      SDF <- SubDF(Exp,site,i,j)
      BDF <- Vienna_BoxDF(Exp, DF=SDF)
      
      plot(SDF$x, SDF$y, type="p", pch=".", col="gray64", xlim=c(0,14), ylim=c(-14,14), cex = 4, 
           xlab="", ylab="", main="", axes=FALSE, yaxs="i")
      axis(side=1, at=seq(0,14,14), cex.axis=1.5, font=2, labels = F)
      axis(side=2, at=seq(-14,14,14), cex.axis=1.5, font=2, labels = F)
      box()
      abline(h=0)
      arrows(x0=BDF$XMean, y0=BDF$YQ1, x1=BDF$XMean, y1=BDF$YQ3, lwd=6, code=0, col="green")
      points(BDF$XMean, BDF$YMedian, pch=16, col="blue", cex=2.5) # Median
      points(BDF$XMean, BDF$YMean, pch=16, col="red", cex=2.5) # Mean
    }
  }
  dev.off()
}

Vienna_SubstitutionTotal_SeqExp("Vienna")

################################################ Vienna_SubstitutionGtoC_Duplex ################################################
Vienna_SubstitutionGtoC_Duplex <- function(){
  SE <- function(Exp,site,i,j){
    if(i == j){return(NA)}
    n <- length(Exp)
    index1 <- which(0:(n-1)%/%4^(9-site)%%4 == i)
    index2 <- which(0:(n-1)%/%4^(9-site)%%4 == j)
    DF <- data.frame(Exp[index1], (Exp[index2]-Exp[index1]))
    colnames(DF) <- c("x","y")
    
    DF <- DF[which(!is.na(DF$x)),]
    DF <- DF[which(!is.na(DF$y)),]
    
    return(DF)
  }
  
  PDF <- function(DF){
    PDF <- data.frame(matrix(NA,100,3))
    colnames(PDF) <- c("XValue","YValue","Counts")
    n <- 0
    for(a in 0:9){
      for(b in -a:(-a+9)){
        n <- n+1
        PDF$XValue[n] <- a
        PDF$YValue[n] <- b
        PDF$Counts[n] <- length(which(DF$x==a & DF$y==b))
      }
    }
    return(PDF)
  }
  
  SDF <- function(DF){
    SDF <- data.frame(matrix(NA,10,5))
    colnames(SDF) <- c("XMean","YMean","YQ1","YMedian","YQ3")
    for(c in 0:9){
      indexBG <- which(DF$x==c)
      SDF$XMean[c+1] <- c
      SDF$YMean[c+1] <- mean(DF$y[indexBG], na.rm = T)
      SDF$YQ1[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.25, type = 2))
      SDF$YMedian[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.5, type=2))
      SDF$YQ3[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.75, type = 2))
    }
    
    return(SDF)
  }
  
  SubstituionViennaLength <- function(f){
    Data <- read.table(file = paste0("../../data/Vienna Package Prediction/",f,"_Duplex.csv"), header = T, sep = ",") # input ordered and full data
    
    i <- 2
    j <- 1
    
    PlotMat <- data.frame(matrix(data = NA, nrow = 10, ncol = 10))
    colnames(PlotMat) <- c("Initial Length",paste0("Length Change (Position",1:9,")"))
    PlotMat$`Initial Length` <- 0:9
    
    for(site in 1:9){
      DF <- SE(Data[,2],site,i,j) # DataFrame
      PDF <- PDF(DF) # Plotting DataFrame
      SDF <- SDF(DF) # Statistic DataFrame
      
      PlotMat[,site+1] <- SDF$YMean
    }
    return(PlotMat)
  }
  
  f = "Vienna"
  output <- SubstituionViennaLength(f)
  write.table(x = output, file = paste0("../../Result/R_Substitution/",f,"_Duplex_total.csv"), quote = F, sep = ",", row.names = F)
}

Vienna_SubstitutionGtoC_Duplex()

################################################ Vienna_SubstitutionTotal_Duplex ################################################
Vienna_SubstitutionTotal_Duplex <- function(){
  SE <- function(Exp,site,i,j){
    if(i == j){return(NA)}
    n <- length(Exp)
    index1 <- which(0:(n-1)%/%4^(9-site)%%4 == i)
    index2 <- which(0:(n-1)%/%4^(9-site)%%4 == j)
    DF <- data.frame(Exp[index1], (Exp[index2]-Exp[index1]))
    colnames(DF) <- c("x","y")
    
    DF <- DF[which(!is.na(DF$x)),]
    DF <- DF[which(!is.na(DF$y)),]
    
    return(DF)
  }
  
  PDF <- function(DF){
    PDF <- data.frame(matrix(NA,100,3))
    colnames(PDF) <- c("XValue","YValue","Counts")
    n <- 0
    for(a in 0:9){
      for(b in -a:(-a+9)){
        n <- n+1
        PDF$XValue[n] <- a
        PDF$YValue[n] <- b
        PDF$Counts[n] <- length(which(DF$x==a & DF$y==b))
      }
    }
    return(PDF)
  }
  
  SDF <- function(DF){
    SDF <- data.frame(matrix(NA,10,5))
    colnames(SDF) <- c("XMean","YMean","YQ1","YMedian","YQ3")
    for(c in 0:9){
      indexBG <- which(DF$x==c)
      SDF$XMean[c+1] <- c
      SDF$YMean[c+1] <- mean(DF$y[indexBG], na.rm = T)
      SDF$YQ1[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.25, type = 2))
      SDF$YMedian[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.5, type = 2))
      SDF$YQ3[c+1] <- as.numeric(quantile(DF$y[indexBG], 0.75, type = 2))
    }
    return(SDF)
  }
  
  SubstitutionEffect <- function(f){
    Data <- read.table(file = paste0("../../data/Vienna Package Prediction/",f,"_Duplex.csv"), header = T, sep = ",") # input ordered and full data
    
    char <- c("A","C","G","U")
    BaseIndex <- matrix(data=c(0,0,0,3,3,3,2,2,2,1,1,1,3,2,1,0,2,1,0,3,1,0,3,2), nrow=12, ncol=2)
    tiff(filename = paste0("../../Result/R_Substitution/",f,"_Duplex_total.tiff"), width = 1800, height = 4800, units = "px")
    par(mai = c(.1, .1, .1, .1), mfrow = c(12,9))
    
    for(k in 1:12){
      i <- BaseIndex[k,1]
      j <- BaseIndex[k,2]
      
      for(site in 1:9){
        DF <- SE(Data[,2],site,i,j) # DataFrame
        PDF <- PDF(DF) # Plotting DataFrame
        SDF <- SDF(DF) # Statistic DataFrame
        
        plot(PDF$XValue, PDF$YValue, pch=20, col="gray64", cex=log10(PDF$Counts)*1.5+0.1, xlim=c(0,9), ylim=c(-9,9), 
             xlab="", ylab="", main="", axes=FALSE)
        axis(side=1, at=0:9, cex.axis=1.5, font=2, labels = F)
        axis(side=2, at=-9:9, cex.axis=1.5, font=2, labels = F)
        box()
        abline(h = 0)
        BarIndex <- which((SDF$YQ3-SDF$YQ1)!=0)
        arrows(x0=SDF$XMean[BarIndex], y0=SDF$YQ1[BarIndex], x1=SDF$XMean[BarIndex], 
               y1=SDF$YQ3[BarIndex], lwd=9, code=0, col="green")
        points(SDF$XMean, SDF$YMedian, pch=16, col="blue", cex=4) # Median
        points(SDF$XMean, SDF$YMean, pch=16, col="red", cex=4) # Mean
      }
    }
    dev.off()
  }
  
  SubstitutionEffect("Vienna")
}

Vienna_SubstitutionTotal_Duplex()


