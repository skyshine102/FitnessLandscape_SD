### ====
NucleotideContent <- function(f,ExpMax=3.2){
  Data <- read.table(file = paste0("../../data/Read Count Table/UnionTable/",f,".csv"), header = T, sep = ",")
  
  CharTable <- data.frame(matrix(data = 0, nrow = 4^9, ncol = 4))
  colnames(CharTable) <- c("A","C","G","U")
  for(i in 0:8){
    CharTable$A[0:(4^9-1) %/% 4^i %% 4 == 0] <- CharTable$A[0:(4^9-1) %/% 4^i %% 4 == 0] + 1
    CharTable$C[0:(4^9-1) %/% 4^i %% 4 == 1] <- CharTable$C[0:(4^9-1) %/% 4^i %% 4 == 1] + 1
    CharTable$G[0:(4^9-1) %/% 4^i %% 4 == 2] <- CharTable$G[0:(4^9-1) %/% 4^i %% 4 == 2] + 1
    CharTable$U[0:(4^9-1) %/% 4^i %% 4 == 3] <- CharTable$U[0:(4^9-1) %/% 4^i %% 4 == 3] + 1
  }
  
  color <- c("#008F00","#4180FF","#000000","#FF2600")
  
  PlotData <- data.frame(matrix(data = NA, nrow = 40, ncol = 3))
  colnames(PlotData) <- c("content","mean","sd")
  PlotData$content <- rep(0:9, 4)
  for(i in 0:3){
    for(j in 1:10){
      PlotData$mean[i*10+j] <- mean(Data$LogMean[CharTable[,i+1] == j-1], na.rm = T)
      PlotData$sd[i*10+j] <- sd(Data$LogMean[CharTable[,i+1] == j-1], na.rm = T)
    }
  }
  PlotData[is.na(PlotData)] <- 0
  
  svg(filename = paste0("../../Result/nucleotide_trend/",f,"_R_NucleotideTrend.svg"), width = 8, height = 6) # open picture file
  
  plot(-10, -10, xlab = "Nucleotide Content", ylab = "expression",
       xlim = c(0,9), ylim = c(0,ExpMax), axes = FALSE, yaxs = "i")
  axis(side = 1, at = c(0:9))
  axis(side = 2, at = seq(0,ExpMax,ExpMax/4))
  box()
  
  for(i in 0:3){
    PlotData$content[1:10+i*10] <- PlotData$content[1:10+i*10]+i*0.1-0.1
    arrows(x0 = PlotData$content[1:9+i*10], y0 = PlotData$mean[1:9+i*10]-PlotData$sd[1:9+i*10], 
           x1 = PlotData$content[1:9+i*10], y1 = PlotData$mean[1:9+i*10]+PlotData$sd[1:9+i*10], 
           length = 0.1, lwd = 2, angle = 90, code = 3, lty = 1, col = color[i+1])
  }
  
  for(i in 0:3){
    lines(x = PlotData$content[1:10+i*10], y = PlotData$mean[1:10+i*10], pch = 20, col = color[i+1], lwd = 5)
  }
  
  legend(1,ExpMax-0.2,c("A","C","G","U"), pch = 20, lty = 1, col = color)
  
  dev.off() # close picture file
}
### ====
NucleotideContent("dmsC_SDR_union_count25")
NucleotideContent("arti_SDR_union_count25")
NucleotideContent("fepB_SDR_union_count25")










