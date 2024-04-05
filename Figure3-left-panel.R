#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the analysis done for Figure-3 left panel.			#
#################################################################################################

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)
source("suppFunctions.R")


#####################################################################--------FIGURE-3 LEFT PANEL-----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure3/Fig3_NBD_LP_table.txt"
otu_t<-read.table(otu_fn, sep=",", header=TRUE, row.names = 1)
otu_t <- filter(otu_t,0.001,0.10)
otu_t


#copy number normalization
otu_t <- rrnNorm (otu_t)
otu_t


#clubbing genera where species are known
otu_t[c("Acutalibacter muris"),]  <- colSums(otu_t[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t[c("Turicimonas muris"),]  <- colSums(otu_t[c("Burkholderiales bacterium YL45"),])
otu_t[c("Clostridium innocuum"),]  <- otu_t[c("Erysipelotrichaceae bacterium I46"),]
otu_t[c("Enterocloster"),]  <- otu_t[c("Lachnoclostridium sp. YL32"),]
otu_t
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis")
otu_t  <- otu_t[!(row.names(otu_t) %in% rem),]

#filter by abundance threshold and prevalence
#otu_t2 <- filter(otu_t,0.001,0.75)

#add unclassified
otu_t2 <- rbind(otu_t2,unclassified=0)
for(i in 1:ncol(otu_t2)) {otu_t2["unclassified",i] <- 1 - sum(otu_t2[,i])}

#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),]) 
stats1
write.csv(stats1, "csv_files/Fig3_NBD_LP_sampleStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig3_NBD_LP")
########################################################################################------------END-------############################################################################################


