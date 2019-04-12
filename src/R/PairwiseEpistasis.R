

############################################## Pairwise epistasis for experimental data ##############################################
### Parameter ====
f <- "dmsC_SDR_union"
spacer <- 4
ExpMax <- 0.14
#f <- "arti_SDR_union"
#spacer <- 5
#ExpMax <- 0.16
f <- "fepB_SDR_union"
spacer <- 6
ExpMax <- 0.34

### Calculator ====
Data <- read.table(file = paste0("../../data/Epistasis Table/",f,"_order2.csv"), header = T, sep = ",")
char <- c("A","U","C","G")
NList <- character()
index <- 0
for(i in (-9-spacer):(-1-spacer)){
  for(j in 1:4){
    index <- index+1
    NList[index] <- paste0(as.character(i),char[j])
  }
}

epi <- matrix(data = 0, nrow = 36, ncol = 36, dimnames = list(NList,NList))
index <- 0
for(sitei in 0:7){
  for(sitej in (sitei+1):8){
    for(ni in 1:4){
      for(nj in 1:4){
        index <- index+1
        epi[sitei*4+ni,sitej*4+nj] <- Data$LogMean[index] # for experimental data
        
      }
    }
  }
}

epi <- epi + t(epi)
epi[epi == 0] <- NA

### Correct Nucleotide Order ====
temp <- epi[seq(2,34,4),]
epi[seq(2,34,4),] <- epi[seq(4,36,4),]
epi[seq(4,36,4),] <- temp

temp <- epi[seq(3,35,4),]
epi[seq(3,35,4),] <- epi[seq(4,36,4),]
epi[seq(4,36,4),] <- temp

temp <- epi[,seq(2,34,4)]
epi[,seq(2,34,4)] <- epi[,seq(4,36,4)]
epi[,seq(4,36,4)] <- temp

temp <- epi[,seq(3,35,4)]
epi[,seq(3,35,4)] <- epi[,seq(4,36,4)]
epi[,seq(4,36,4)] <- temp

### ====
NNpattern <- matrix(data = 0, nrow = 4, ncol = 4, dimnames = list(char,char))
for(i in seq(1,36,4)){
  if(i+7 < 36){
    NNpattern <- NNpattern + epi[i:(i+3),(i+4):(i+7)]
  }
  if(i+11 < 36){
    NNpattern <- NNpattern + epi[i:(i+3),(i+8):(i+11)]
  }
}
NNpattern <- NNpattern + t(NNpattern)
NNpattern <- NNpattern/30 * 1.3 # correct the color scale

### Plot ====
library(lattice)

svg(filename = paste0("../../Result/PairwiseEpistasis/",f,"_36by36.svg"), width = 8, height = 6) # open picture file
lattice.options(axis.padding=list(factor=0.5))
rgb.palette <- colorRampPalette(c("blue","white", "red"), space = "rgb")
levelplot(epi, main="", xlab="", ylab="", 
          col.regions=rgb.palette(100), cuts=100, at=c(seq(-ExpMax, ExpMax, length.out = 101)), 
          scales=list(y=(list()), tck = c(1,0), x=list(rot=90)),
          colorkey=list(labels=list(at=seq(-ExpMax, ExpMax, length.out = 5))),
          panel = function(...){
            panel.fill(col = "gray")
            panel.levelplot(...)
            panel.abline(h = seq(4.5,32.5,4), v = seq(4.5,32.5,4))
          })
dev.off() # close picture file

svg(filename = paste0("../../Result/PairwiseEpistasis/",f,"_4by4.svg"), width = 8, height = 6) # open picture file
lattice.options(axis.padding=list(factor=0.5))
rgb.palette <- colorRampPalette(c("blue","white", "red"), space = "rgb")
levelplot(NNpattern, main="", xlab="", ylab="", 
          col.regions=rgb.palette(100), cuts=100, 
          scales=list(y=(list()), tck = c(1,0), x=list(rot=0)), 
          colorkey=list(labels=list(at=seq(-ExpMax, ExpMax, length.out = 5))), at=seq(-ExpMax, ExpMax, length.out = 101), 
          panel = function(...){
            panel.fill(col = "gray")
            panel.levelplot(...)
          })
dev.off() # close picture file

############################################## Pairwise epistasis for Vienna data ##############################################
### Parameter ====
f <- "Vienna"
spacer <- -10
ExpMax <- 1

### Calculator ====
Data <- read.table(file = paste0("../../data/Vienna Package Prediction/",f,"_order2.csv"), header = T, sep = ",")
char <- c("A","U","C","G")
NList <- character()
index <- 0
for(i in (-9-spacer):(-1-spacer)){
  for(j in 1:4){
    index <- index+1
    NList[index] <- paste0(as.character(i),char[j])
  }
}

epi <- matrix(data = 0, nrow = 36, ncol = 36, dimnames = list(NList,NList))
index <- 0
for(sitei in 0:7){
  for(sitej in (sitei+1):8){
    
    for(ni in 1:4){
      for(nj in 1:4){
        index <- index+1
        epi[sitei*4+ni,sitej*4+nj] <- Data$LogMean[index] # for Vienna data
        
      }
    }
  }
}

epi <- epi + t(epi)
epi[epi == 0] <- NA

### Correct Nucleotide Order ====
temp <- epi[seq(2,34,4),]
epi[seq(2,34,4),] <- epi[seq(4,36,4),]
epi[seq(4,36,4),] <- temp

temp <- epi[seq(3,35,4),]
epi[seq(3,35,4),] <- epi[seq(4,36,4),]
epi[seq(4,36,4),] <- temp

temp <- epi[,seq(2,34,4)]
epi[,seq(2,34,4)] <- epi[,seq(4,36,4)]
epi[,seq(4,36,4)] <- temp

temp <- epi[,seq(3,35,4)]
epi[,seq(3,35,4)] <- epi[,seq(4,36,4)]
epi[,seq(4,36,4)] <- temp

### ====
NNpattern <- matrix(data = 0, nrow = 4, ncol = 4, dimnames = list(char,char))
for(i in seq(1,36,4)){
  if(i+7 < 36){
    NNpattern <- NNpattern + epi[i:(i+3),(i+4):(i+7)]
  }
  if(i+11 < 36){
    NNpattern <- NNpattern + epi[i:(i+3),(i+8):(i+11)]
  }
}
NNpattern <- NNpattern + t(NNpattern)
NNpattern <- NNpattern/30

### Plot ====
library(lattice)

svg(filename = paste0("../../Result/PairwiseEpistasis/",f,"_36by36.svg"), width = 8, height = 6) # open picture file
lattice.options(axis.padding=list(factor=0.5))
rgb.palette <- colorRampPalette(c("red","white", "blue"), space = "rgb") # Vienna
levelplot(epi, main="", xlab="", ylab="", 
          col.regions=rgb.palette(100), cuts=100, at=c(seq(-ExpMax, ExpMax, length.out = 101)), 
          scales=list(y=(list()), tck = c(1,0), x=list(rot=90)),
          colorkey=list(labels=list(at=seq(-ExpMax, ExpMax, length.out = 5))),
          panel = function(...){
            panel.fill(col = "gray")
            panel.levelplot(...)
            panel.abline(h = seq(4.5,32.5,4), v = seq(4.5,32.5,4))
          })
dev.off() # close picture file

svg(filename = paste0("../../Result/PairwiseEpistasis/",f,"_4by4.svg"), width = 8, height = 6) # open picture file
lattice.options(axis.padding=list(factor=0.5))
rgb.palette <- colorRampPalette(c("red","white", "blue"), space = "rgb")
levelplot(NNpattern, main="", xlab="", ylab="", 
          col.regions=rgb.palette(100), cuts=100, 
          scales=list(y=(list()), tck = c(1,0), x=list(rot=0)), 
          colorkey=list(labels=list(at=seq(-ExpMax, ExpMax, length.out = 5))), at=seq(-ExpMax, ExpMax, length.out = 101), 
          panel = function(...){
            panel.fill(col = "gray")
            panel.levelplot(...)
          })
dev.off() # close picture file
