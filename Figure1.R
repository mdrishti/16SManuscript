#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the analysis done for Figure-1.				#
#################################################################################################

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)
source("suppFunctions.R")


#####################################################################--------FIGURE-1 LEFT PANEL -----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure1/Fig1_LP_table.txt"
otu_t <- read.table(otu_fn, sep=",", header=TRUE, row.names = 1)
otu_t


#copy number normalization
otu_t2 <- rrnNorm (otu_t)
otu_t2


#clubbing genera where species are known
#otu_t2[c("Acutalibacter muris"),]  <- colSums(otu_t2[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t2[c("Turicimonas muris"),]  <- colSums(otu_t2[c("Burkholderiales bacterium YL45"),])
otu_t2[c("Clostridium innocuum"),]  <- otu_t2[c("Erysipelotrichaceae bacterium I46"),]
otu_t2[c("Enterocloster"),]  <- colSums(otu_t2[c("Lachnoclostridium sp. YL32"),])
otu_t2
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis")
otu_t2  <- otu_t2[!(row.names(otu_t2) %in% rem),]
otu_t2


#filter by abundance threshold and prevalence
#otu_t2 <- filter(otu_t2,0.001,0.75) ###########################
#otu_t2 #########################################################


#add unclassified
otu_t2 <- rbind(otu_t2,unclassified=0)
for(i in 1:ncol(otu_t2)) {otu_t2["unclassified",i] <- 1 - sum(otu_t2[,i])}


#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),]) 
stats1
write.csv(stats1, "csv_files/Fig1_LP_sampleStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig1_LP")
########################################################################################------------END-------############################################################################################




#####################################################################--------FIGURE-1 RIGHT PANEL -----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure1/Fig1_RP_table.txt"
otu_t<-read.table(otu_fn, sep=",", header=TRUE, row.names = 1)
otu_t <- filter(otu_t,0,0.75)
otu_t


#copy number normalization
#otu_t2 <- rrnNorm (otu_t)
otu_t2  <- otu_t #No copy number normalization because it is not done for the Illumina data


#clubbing genera where species are known
otu_t2[c("Acutalibacter muris"),]  <- colSums(otu_t2[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t2[c("Turicimonas muris"),]  <- colSums(otu_t2[c("Burkholderiales bacterium YL45"),])
otu_t2[c("Clostridium innocuum"),]  <- otu_t2[c("Erysipelotrichaceae bacterium I46"),]
otu_t2[c("Enterocloster"),]  <- colSums(otu_t2[c("Lachnoclostridium sp. YL32"),])
#otu_t2[c("Blautia"),]  <- colSums(otu_t2[c("Blautia sp. YL58","Blautia producta","Blautia hominis"),])
otu_t2
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis","Total_NF")
otu_t2  <- otu_t2[!(row.names(otu_t2) %in% rem),]
otu_t2


#filter by abundance threshold and prevalence
otu_t2 <- filter(otu_t2,0.001,0.10)
otu_t2


#add unclassified
otu_t2 <- rbind(otu_t2,unclassified=0)
for(i in 1:ncol(otu_t2)) {otu_t2["unclassified",i] <- 1 - sum(otu_t2[,i])}
otu_t2

#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),]) 
stats1
write.csv(stats1, "csv_files/Fig1_RP_sampleStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig1_RP")
########################################################################################------------END-------############################################################################################


