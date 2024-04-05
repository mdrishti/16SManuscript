#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the analysis done for Figure-2.				#
#################################################################################################

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)
source("suppFunctions.R")


#####################################################################--------FIGURE-2 RIGHT PANEL (APPROACH-2)-----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure2/Fig2_approach2_RP_table.txt"
otu_t<-read.table(otu_fn, sep="\t", header=TRUE, row.names = 1)


#copy number normalization
otu_t2 <- rrnNorm (otu_t)

#clubbing genera where species are known
otu_t2[c("Acutalibacter muris"),]  <- colSums(otu_t2[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t2[c("Turicimonas muris"),]  <- colSums(otu_t2[c("Burkholderiales bacterium YL45"),])
otu_t2[c("Clostridium innocuum"),]  <- otu_t2[c("Erysipelotrichaceae bacterium I46"),]
otu_t2[c("Enterocloster"),]  <- colSums(otu_t2[c("Lachnoclostridium sp. YL32"),])
otu_t2
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis")
otu_t2  <- otu_t2[!(row.names(otu_t2) %in% rem),]
otu_t2

#filter by abundance threshold and prevalence
otu_t2 <- filter(otu_t2,0.001,0.10)
otu_t2

#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)[which(!rownames(otu_t2) %in% c("[Clostridium] scindens"))]
expecT2  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),c("barcode01","barcode02","barcode03","barcode04","barcode05")]) 
stats2 <- preRec(expecT2, otu_t2[which(rownames(otu_t2) != "unclassified"),c("barcode06","barcode07","barcode09","barcode10")]) 
stats1
stats2
write.csv(rbind(stats1,stats2), "csv_files/Fig2_NBD_RP_sampleStats.csv")


#calculate L1-norm, L2-norm, precision and recall for mock
mockD <- readMock()
mock <- c("Bacteroides caecimuris","Bifidobacterium animalis","Blautia sp. YL58","Turicimonas muris","Enterococcus faecalis","Clostridium innocuum","Flavonifractor plautii","Enterocloster","Lactobacillus reuteri")
b11 <- otu_t2[rownames(otu_t2) %in%mock,"barcode11", drop=FALSE]
b12 <- otu_t2[rownames(otu_t2) %in%mock,"barcode12", drop=FALSE]
b13 <- otu_t2[rownames(otu_t2) %in%mock,"barcode13", drop=FALSE]
m1 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b11[match(mock, rownames(b11)),,drop=FALSE])
m2 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b12[match(mock, rownames(b12)),,drop=FALSE])
m3 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b13[match(mock, rownames(b13)),,drop=FALSE])
m <- data.frame(m1=m1,m2=m2,m3=m3)
m
mock <- c("Bacteroides caecimuris","Bifidobacterium animalis","Blautia sp. YL58","Turicimonas muris","Enterococcus faecalis","Clostridium innocuum","Flavonifractor plautii","Enterocloster","Lactobacillus reuteri", "L1-norm", "L2-norm", "True positives", "precision", "recall")
rownames(m) <- mock
write.csv(m,"csv_files/Fig2_NBD_RP_mockStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig2_NBD_RP")
########################################################################################------------END-------############################################################################################



#####################################################################--------FIGURE-2 LEFT PANEL (APPROACH-1)-----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure2/Fig2_approach1_LP_table.txt"
otu_t<-read.table(otu_fn, sep="\t", header=TRUE, row.names = 1)
otu_Bif <- otu_t["Bifidobacterium animalis",] #for later, the relative abundance is zero and thus gets filtered in the next steps


#copy number normalization
otu_t2 <- rrnNorm (otu_t)

#clubbing genera where species are known
otu_t2[c("Acutalibacter muris"),]  <- colSums(otu_t2[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t2[c("Turicimonas muris"),]  <- colSums(otu_t2[c("Burkholderiales bacterium YL45"),])
otu_t2[c("Clostridium innocuum"),]  <- otu_t2[c("Erysipelotrichaceae bacterium I46"),]
otu_t2[c("Enterocloster"),]  <- colSums(otu_t2[c("Lachnoclostridium sp. YL32"),])
otu_t2
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis")
otu_t2  <- otu_t2[!(row.names(otu_t2) %in% rem),]
otu_t2


#filter by abundance threshold and prevalence
otu_t2 <- filter(otu_t2,0.001,0.25)
otu_t2 <- rbind(otu_t2,otu_Bif)


#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)[which(!rownames(otu_t2) %in% c("[Clostridium] scindens"))]
expecT2  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),c("barcode01","barcode02","barcode04","barcode05","barcode06")]) 
stats2 <- preRec(expecT2, otu_t2[which(rownames(otu_t2) != "unclassified"),c("barcode12","barcode13","barcode15","barcode16")]) 
stats1
stats2
write.csv(rbind(stats1,stats2), "csv_files/Fig2_16SSQK024_LP_sampleStats.csv")


mock <- c("Bacteroides caecimuris","Bifidobacterium animalis","Blautia sp. YL58","Turicimonas muris","Enterococcus faecalis","Clostridium innocuum","Flavonifractor plautii","Enterocloster","Lactobacillus reuteri")
b21 <- otu_t2[rownames(otu_t2) %in%mock,"barcode21", drop=FALSE]
b22 <- otu_t2[rownames(otu_t2) %in%mock,"barcode22", drop=FALSE]
b23 <- otu_t2[rownames(otu_t2) %in%mock,"barcode23", drop=FALSE]
m1 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b21[match(mock, rownames(b21)),,drop=FALSE])
m2 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b22[match(mock, rownames(b22)),,drop=FALSE])
m3 <- eucl.vec(mockD[,"rrn_mod",drop=FALSE],b23[match(mock, rownames(b23)),,drop=FALSE])
m <- data.frame(m1=m1,m2=m2,m3=m3)
m
mock <- c("Bacteroides caecimuris","Bifidobacterium animalis","Blautia sp. YL58","Turicimonas muris","Enterococcus faecalis","Clostridium innocuum","Flavonifractor plautii","Enterocloster","Lactobacillus reuteri", "L1-norm", "L2-norm", "True positives", "precision", "recall")
rownames(m) <- mock
write.csv(m,"csv_files/Fig2_16SSQK024_LP_mockStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig2_16SSQK024_LP")
########################################################################################------------END-------############################################################################################
