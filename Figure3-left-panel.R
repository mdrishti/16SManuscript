#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the analysis done for Figure-3 left panel.				#
#################################################################################################

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)
source("suppFunctions.R")


#####################################################################--------FIGURE-3 LEFT PANEL-----------##############################################################
#load the OTU table and the count table:
otu_fn <- "data/Figure3/NBD_LeftPanel_Table_all"
otu_t<-read.table(otu_fn, sep="\t", header=TRUE, row.names = 1)
#for(i in 1:ncol(otu_t)){print(sum(otu_t[1:(nrow(otu_t)),i]))} #check
#length(which(rowSums(otu_t) != 0)) #check


#for barcode fitering----remove after extracting table
bc <- c("barcode49","barcode50","barcode51","barcode52","barcode53")
otu_t <- otu_t[,colnames(otu_t) %in% bc]


#copy number normalization
otu_t2 <- rrnNorm (otu_t)

#clubbing genera where species are known
otu_t[c("Acutalibacter muris"),]  <- colSums(otu_t[c("Acutalibacter muris","Hungateiclostridiaceae bacterium KB18"),])
otu_t[c("Turicimonas muris"),]  <- colSums(otu_t[c("Burkholderiales bacterium YL45"),])
otu_t[c("Clostridium innocuum"),]  <- otu_t[c("Erysipelotrichaceae bacterium I46"),]
otu_t[c("Enterocloster"),]  <- colSums(otu_t[c("Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis"),])
otu_t
rem  <- c("Hungateiclostridiaceae bacterium KB18","Burkholderiales bacterium YL45","Turicimonas muris YL45","Erysipelotrichaceae bacterium I46","Lachnoclostridium sp. YL32","Enterocloster bolteae","Enterocloster clostridioformis")
otu_t  <- otu_t[!(row.names(otu_t) %in% rem),]
otu_t

#filter by abundance threshold and prevalence
otu_t2 <- filter(otu_t,0.001,0.75)
otu_t2

#add unclassified
otu_t2 <- rbind(otu_t2,unclassified=0)
for(i in 1:ncol(otu_t2)) {otu_t2["unclassified",i] <- 1 - sum(otu_t2[,i])}

#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2)   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),]) 
stats1
write.csv(stats1, "NBD_Figure3_Samples_stats_X.csv")


#Plot barplots
meltAndPlot(otu_t2, "NBD_Fig3")
########################################################################################------------END-------############################################################################################


