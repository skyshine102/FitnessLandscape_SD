### SingleNucleotideEffect ====
SingleNucleotideEffect <- function(f,spacing = 5, ExpMax = 0.34){
  Data <- read.table(file = paste0("../../data/Epistasis Table/",f,"_order1.csv"), sep = ",", header = T)
  
  svg(filename = paste0("SingleNucleotide/",f,"_SingleNucleotide.svg"), width = 8, height = 6) # open file
  
  color <- c("#008F00","#4180FF","#000000","#FF2600")
  plot(0,0,type='n', xlim = c(-16,-4), ylim = c(-ExpMax,ExpMax), axes = F, xaxs = "i", yaxs = "i", xlab = "Position", ylab = "Effect", main = "") # for experimental data
  #plot(0,0,type='n', xlim = c(0,10), ylim = c(-ExpMax,ExpMax), axes = F, xaxs = "i", yaxs = "i", xlab = "Position", ylab = "Effect", main = "") # for Vienna data
  axis(side = 1, at = seq(-16,-4,2)) # for experimental data
  #axis(side = 1, at = seq(0,10,2)) # for Vienna data
  axis(side = 2, at = seq(-ExpMax,ExpMax,ExpMax/2))
  box()
  for(i in 0:3){
    lines(-(spacing+9):-(spacing+1), Data$LogMean[0:35 %% 4 == i], type = "o", pch = 20, col = color[i+1], lwd = 3) # for experimental data
    #lines(1:9, -Data$dG[0:35 %% 4 == i], type = "o", pch = 20, col = color[i+1], lwd = 3) # for Vienna data
  }
  legend(-5.5, ExpMax-0.01, legend = c("A","C","G","U"), col = color, pch = 20) # for experimental data
  #legend(8.5, -0.4, legend = c("A","C","G","U"), col = color, pch = 20) # for Vienna data
  
  dev.off() # close file
}

### ====
SingleNucleotideEffect("dmsC_SDR_union", spacing = 4)
SingleNucleotideEffect("arti_SDR_union", spacing = 5)
SingleNucleotideEffect("fepB_SDR_union", spacing = 6)
#SingleNucleotideEffect("Vienna", ExpMax = 1.2)
