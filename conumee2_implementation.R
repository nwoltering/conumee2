# for loading 450/850k
library(minfi)
#for loading 935k
library(illuminaio)
library(tidyverse)
# plotly has to be installed before installing conumee2
library(plotly)
library(conumee2.0)
library(conumee)
library(openxlsx) 
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

#450/850k Loading------------------------

#------------------------------
# Sample-Sheet und Rohdaten einlesen
#------------------------------


#----------------------------Einlesen XLSX
samples <- as.data.frame(read.xlsx("C:\\Daten\\HMB\\idat\\SampleSheet.xlsx",rowNames = T))
#samples <- as.data.frame(read.xlsx("/home/niklas/Daten/HMB_MS/SampleSheet_Klinisch.xlsx",rowNames = T))
#samples <- tibble::rownames_to_column(samples, "Sentrix_ID")
samples <- cbind.data.frame(samples, Basename=paste0("C:\\Daten\\HMB\\idat\\",samples$Sentrix_ID))
#samples <- cbind.data.frame(samples, Basename=paste0("/home/niklas/Daten/HMB_MS/",samples$Sentrix_ID))

mit450k <- FALSE


if(mit450k == FALSE){
  RGset<- read.metharray.exp(targets=targets, verbose=TRUE, force = TRUE)
  
  
}else{
  targets_450K <- samples[samples$Panel=="IlluminaHumanMethylation450k",]
  targets_EPIC <- samples[samples$Panel=="IlluminaHumanMethylationEPIC",]
  targets <- rbind(targets_450K,targets_EPIC)
  # Kombinieren der unterschiedlichen Analyseformate
  RGset_450K <- read.metharray.exp(targets=targets_450K, force = TRUE)
  RGset_EPIC <- read.metharray.exp(targets=targets_EPIC, verbose=TRUE, force = TRUE)
  RGset <- combineArrays(RGset_450K,RGset_EPIC) #Zusammenfügen der Arrays, nach Probe-Namen, die Funktion unterstützt max. 2
  rm(RGset_450K,RGset_EPIC)
}
#options(encoding = "UTF-8") 
sampleNames(RGset) <- targets$Name

detP <- detectionP(RGset)


# remove poor quality samples
keep <- colMeans(detP) < 0.05
RGset <- RGset[,keep]

detP <- detP[,keep]


#mSetSq <- preprocessQuantile(RGset)#PreprocessQuantile kann manchmal komische Sachen machen 
mSetSq <- preprocessIllumina(RGset)
mSetSq <- mapToGenome(mSetSq)
mSetSq <- ratioConvert(mSetSq)
mSetRaw <- preprocessRaw(RGset)


# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(RGset, sampGroups=targets$VHL,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$VHL)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$VHL,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$VHL)), 
       text.col=brewer.pal(8,"Dark2"))

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]


# if your data includes males and females, remove probes on the sex chromosomes
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
#                                                       c("chrX","chrY")])
# table(keep)
# mSetSqFlt <- mSetSqFlt[keep,]


# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)

bVals <- getBeta(mSetSqFlt)



GMset <- mapToGenome(mSetRaw)

dropSexChr <- function(object){
  xyIndex <- c(which(seqnames(object) == "chrX"),which(seqnames(object) == "chrY"))
  object[-xyIndex,]
}

GMset <- dropSexChr(GMset)
#Es müssen mal wieder die Chromosomen gefiltert werden

GMset <- addSnpInfo(GMset) #Annotation und entfernen der SNPs
GMset <- dropLociWithSnps(GMset, snps=c("SBE","CpG"), maf=0)
#---------------------------------
# recalculate betas, illumina like
#---------------------------------
methy <- minfi::getMeth(GMset)
unmethy <- getUnmeth(GMset)
betas <- methy / (methy +unmethy +100)
#935k-Loading------------------------
#TODO
#Loading flat reference-------------------
load("C:\\Daten\\HMB\\CNanalysis5_conumee_REF.2017-02-10.RData") 
#Conumee-----------------


#Visualization----------------
