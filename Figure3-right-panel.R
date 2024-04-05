#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the analysis done for Figure-3 right panel.			#
#################################################################################################

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)
source("suppFunctions.R")

#####################################################################--------FIGURE-3 RIGHT PANEL-----------##############################################################
#load the OTU table, taxonomic metadata and sample metdata tables
otu_fn <- "data/Figure3/Fig3_Illumina_RP_table.txt"
otu_t<-read.table(otu_fn, sep=",", header=TRUE, row.names = 1)

#normalization by count
otu_t <- apply(otu_t, 2, function(x) (x/tail(x,1)))

#filter by abundance threshold and prevalence
otu_t2 <- filter(otu_t,0.001,0.75)
otu_t2

#add unclassified
otu_t2 <- rbind(otu_t2,unclassified=0)
for(i in 1:ncol(otu_t2)) {otu_t2["unclassified",i] <- otu_t2["Total",i] - sum(otu_t2[1:(nrow(otu_t2)-2),i])}
otu_t2

#Remove the row labelled "Total"
rem  <- c("Total")
otu_t2  <- data.frame(otu_t2[!(row.names(otu_t2) %in% rem),])
otu_t2

#calculate precision and recall for samples without mock
expecT1  <- rownames(otu_t2[which(rownames(otu_t2) != "unclassified"),])   
stats1 <- preRec(expecT1, otu_t2[which(rownames(otu_t2) != "unclassified"),]) 
stats1
write.csv(stats1, "csv_files/Fig3_Illumina_RP_sampleStats.csv")


#Plot barplots
meltAndPlot(otu_t2, "Figures/Fig3_Illumina_RP")
##################################################################---END---######################################################
