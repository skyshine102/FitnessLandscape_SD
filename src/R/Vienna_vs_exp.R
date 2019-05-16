################################################ dG vs. experimental data ################################################ 
### ====
ViennaCompare <- function(f){
  cex = 3
  
  Model <- read.table(file = "../../data/Vienna Package Prediction/Vienna_SeqExp.csv", header = T, sep = ",")
  Experimental <- read.table(file = paste0("../../data/Read Count Table/UnionTable/",f,".csv"), header = T, sep = ",")
  
  x <- -Model$LogMean
  y <- Experimental$LogMean
  
  min(x)
  max(x)
  min(y, na.rm = T)
  max(y, na.rm = T)
  
  n <- 20
  BoxDF <- data.frame(matrix(data = NA, nrow = n, ncol = 6))
  colnames(BoxDF) <- c("XMean","YMean", "YSD","YQ1","YMedian","YQ3")
  BinList <- seq(0,14,14/n)
  for(i in 1:n){
    index <- which(BinList[i]<x & x<BinList[i+1])
    BoxDF$XMean[i] <- mean(x[index], na.rm = T)
    BoxDF$YMean[i] <- mean(y[index], na.rm = T)
    BoxDF$YSD[i] <- sd(y[index], na.rm = T)
    BoxDF$YQ1[i] <- as.numeric(quantile(y[index], 0.25, type = 2, na.rm = T))
    BoxDF$YMedian[i] <- as.numeric(quantile(y[index], 0.5, type = 2, na.rm = T))
    BoxDF$YQ3[i] <- as.numeric(quantile(y[index], 0.75, type = 2, na.rm = T))
  }
  
  tiff(filename = paste0("../../Result/R_Vienna_vs_exp/",f,"_dG_vs_exp.tiff"), width = 360, height = 360, units = "px")
  
  plot(x, y, pch=".", xlab="", ylab="", xlim=c(0,14), ylim=c(0,3.2), axes=FALSE, yaxs="i", col="gray64", cex=cex)
  axis(side = 1, at = seq(0,14,14/5), cex.axis = 1.5)
  axis(side = 2, at = seq(0,3.2,3.2/4), cex.axis = 1.5)
  box()
  arrows(x0 = BoxDF$XMean, y0 = BoxDF$YMean-BoxDF$YSD, x1 = BoxDF$XMean, BoxDF$YMean+BoxDF$YSD, 
         length = .06, lwd = 2, angle = 90, code = 3, col = "red") # Mean
  points(x = BoxDF$XMean, y = BoxDF$YMean, pch = 20, col = "red", cex = 2.5) # Mean
  
  print(cor.test(x, y, method = "pearson", exact = FALSE)$estimate)
  print(cor.test(x, y, method = "pearson", exact = FALSE)$p.value)
  print(cor.test(x, y, method = "spearman", exact = FALSE)$estimate)
  print(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  
  write.table(x = data.frame(x,y), file = paste0("../../Result/R_Vienna_vs_exp/",f,"_dG_vs_exp.csv"), sep = ",", row.names = F)
  dev.off()
}

### ====
ViennaCompare("arti_SDR_union_count25")
ViennaCompare("dmsC_SDR_union_count25")
ViennaCompare("fepB_SDR_union_count25")

################################################ duplex length vs. experimental data ################################################
### choose bin size and get the information of each bin ====
BinDataFrame <- function(x, y, xlim, n){
  BDF <- data.frame(matrix(data = NA, nrow = n, ncol = 6))
  colnames(BDF) <- c("XMean","YMean", "YSD","YQ1","YMedian","YQ3")
  BinList <- seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/n)
  for(i in 1:n){
    index <- which(BinList[i]<=x & x<BinList[i+1])
    BDF$XMean[i] <- mean(x[index], na.rm = T)
    BDF$YMean[i] <- mean(y[index], na.rm = T)
    BDF$YSD[i] <- sd(y[index], na.rm = T)
    BDF$YQ1[i] <- as.numeric(quantile(y[index], 0.25, type = 2, na.rm = T))
    BDF$YMedian[i] <- as.numeric(quantile(y[index], 0.5, type = 2, na.rm = T))
    BDF$YQ3[i] <- as.numeric(quantile(y[index], 0.75, type = 2, na.rm = T))
  }
  return(BDF)
}

### main function ====
## data processing
ViennaDuplex <- function(f){
  Exp <- read.table(file = paste0("../../data/Read Count Table/UnionTable/",f,".csv"), header = T, sep = ",")
  duplex <- read.table(file = "../../data/Vienna Package Prediction/Vienna_Duplex.csv", header = T, sep = ",")
  
  x <- Exp$LogMean
  y <- duplex[,2]
  y <- y[!is.na(x)]
  x <- x[!is.na(x)]
  XMin <- 3.2
  BDF <- BinDataFrame(x = x, y = y, xlim = c(0,XMin), n = 20)
  
  print(cor.test(x, y, method = "pearson", exact = FALSE)$estimate)
  print(cor.test(x, y, method = "pearson", exact = FALSE)$p.value)
  write.table(x = data.frame(x,y), file = paste0("../../Result/R_Vienna_vs_exp/",f,"_duplex_vs_exp.csv"), sep = ",", row.names = F)
  
  ## plotting
  svg(filename = paste0("../../Result/R_Vienna_vs_exp/",f,"_duplex_vs_exp.svg"), width = 6, height = 6) # open picture file
  
  plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,XMin), ylim=c(-0.5,9.5), axes=FALSE)
  axis(side=1, at=seq(0,XMin,XMin/4), cex.axis=1.5)
  axis(side=2, at=seq(0,9,1), cex.axis=1.5)
  box()
  
  library(vioplot)
  vioplot(x[y==0], col="#e3e4e5", horizontal=TRUE, at=0, add=TRUE, lwd=1, lty=1, lineCol="gray", border="gray")
  for(i in 2:9){
    vioplot(x[y==i], col="#e3e4e5", horizontal=TRUE, at=i, add=TRUE, lwd=1, lty=1, lineCol="gray", border="gray")}
  
  arrows(x0=BDF$XMean, y0=BDF$YMean-BDF$YSD, x1=BDF$XMean, BDF$YMean+BDF$YSD, length=.06, lwd=2, angle=90, code=3, col="red")
  points(x=BDF$XMean, y=BDF$YMean, pch=20, col="red", cex=2.5)
  
  dev.off() # close picture file
}

### ====
ViennaDuplex("arti_SDR_union_count25")
ViennaDuplex("dmsC_SDR_union_count25")
ViennaDuplex("fepB_SDR_union_count25")

################################################ dG vs. duplex length ################################################
### choose bin size and get the information of each bin ====
BinDataFrame <- function(x, y, xlim, n){
  BDF <- data.frame(matrix(data = NA, nrow = n, ncol = 6))
  colnames(BDF) <- c("XMean","YMean", "YSD","YQ1","YMedian","YQ3")
  BinList <- seq(xlim[1], xlim[2], (xlim[2]-xlim[1])/n)
  for(i in 1:n){
    index <- which(BinList[i]<=x & x<BinList[i+1])
    BDF$XMean[i] <- mean(x[index], na.rm = T)
    BDF$YMean[i] <- mean(y[index], na.rm = T)
    BDF$YSD[i] <- sd(y[index], na.rm = T)
    BDF$YQ1[i] <- as.numeric(quantile(y[index], 0.25, type = 2, na.rm = T))
    BDF$YMedian[i] <- as.numeric(quantile(y[index], 0.5, type = 2, na.rm = T))
    BDF$YQ3[i] <- as.numeric(quantile(y[index], 0.75, type = 2, na.rm = T))
  }
  return(BDF)
}

### main function ====
## data processing
dG <- read.table(file = "../../data/Vienna Package Prediction/Vienna_SeqExp.csv", header = T, sep = ",")
duplex <- read.table(file = "../../data/Vienna Package Prediction/Vienna_Duplex.csv", header = T, sep = ",")

x <- -dG$LogMean
y <- duplex[,2]
BDF <- BinDataFrame(x = x, y = y, xlim = c(0,14), n = 20)

print(cor.test(x, y, method = "pearson", exact = FALSE)$estimate)
print(cor.test(x, y, method = "pearson", exact = FALSE)$p.value)
write.table(x = data.frame(x,y), file = paste0("../../Result/R_Vienna_vs_exp/dG_vs_duplex.csv"), sep = ",", row.names = F)

## plotting
svg(filename = "../../Result/R_Vienna_vs_exp/dG_vs_duplex.svg", width = 6, height = 6) # open picture file

plot(0, 0, type="n", xlab="", ylab="", xlim=c(0,14), ylim=c(-0.5,9.5), axes=FALSE)
axis(side=1, at=seq(0,14,14/5), cex.axis=1.5)
axis(side=2, at=seq(0,9,1), cex.axis=1.5)
box()

library(vioplot)
vioplot(c(x[y==0],0.001), col="gray", horizontal=TRUE, at=0, add=TRUE, lwd=5, border="gray")
for(i in 2:9){vioplot(x[y==i], col="#e3e4e5", horizontal=TRUE, at=i, add=TRUE, lwd=1, lty=1, lineCol="gray", border="gray")}

arrows(x0=BDF$XMean, y0=BDF$YMean-BDF$YSD, x1=BDF$XMean, BDF$YMean+BDF$YSD, length=.06, lwd=2, angle=90, code=3, col="red")
points(x=BDF$XMean, y=BDF$YMean, pch=20, col="red", cex=2.5)

dev.off() # close picture file












