#################################################################################################
# @Project - 16S rRNA gene amplicon nanopore sequencing analysis.				#
# @Description - This file gives the basic functions to be used.				#
#################################################################################################


library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(Polychrome)

### copy number normalization
rrnNorm <- function(otu_t) {
	rrn.meta <- "/home/drishti/Softwares/Manuscripts/16S_manuscript/rrn_metadata"
	rrn.metad <- read.table(rrn.meta, sep="\t", header=TRUE, row.names = 1)
	otu_t1 <- apply(otu_t, 2, function(x) (x*tail(x,1))) #Multiply and get actual abundance
	otu_t1 <- data.frame(otu_t1[-c(nrow(otu_t1)),])
	rrn.x <- unlist(rrn.metad[rownames(otu_t1),"rrn"])
	otu_t1$rrn <- rrn.x
	otu_t2 <- apply(otu_t1, 1, function(x) (x/tail(x,1))) #Normalize by copy number
	otu_t2
	otu_t2 <- data.frame(otu_t2[-c(nrow(otu_t2)),])
	otu_t2
	otu_t2[,"Total_NF"]  <-  rowSums(otu_t2)
	otu_t2
	otu_t2 <- apply(otu_t2, 1, function(x) (x/tail(x,1))) #Divide and get relative abundance
	otu_t2
	otu_t2 <- data.frame(otu_t2[-c(nrow(otu_t2)),])
	otu_t2
	rownames(otu_t2) <- rownames(otu_t1)
	otu_t2
	return(otu_t2)
}

#filtering by the presence of the species in atleast x% samples
filter <- function(otu_t1,abThresh,numS) {
	otu_t2 <- otu_t1[rowSums(otu_t1 > abThresh) >= (numS*ncol(otu_t1)),] #Filtering at 75%
	return(otu_t2)
}

#calculate precision and recall for biological samples
preRec <- function(x1,x2) {
        tp <- c()                                                     
        recall <- c()
        precision <- c()
        for(i in 1:ncol(x2)) {tp <- append(tp,length(which( x2[,i] >= 0.01 & rownames(x2[,i, drop=FALSE]) %in% x1== TRUE)))}
        for(i in 1:ncol(x2)) {recall <- append(recall,(tp[i]/length(x1)))}
        for(i in 1:ncol(x2)) {precision <- append(precision,(tp[i]/length(which(x2[,i] != 0))))}
        list1 <- data.frame(tp=tp, precision=precision, recall=recall)
        return(list1)
}

#calculate L1-norm, L2-norm, precision, reacall for mock community
meltAndPlot <- function(otu_t2, prefix) {
	otu_t2$taxa <- rownames(otu_t2)
	otu_t.m  <-  melt(otu_t2, id=c("taxa"))
	print(otu_t.m)
	i <- sapply(otu_t.m, is.factor)
	otu_t.m[i] <- lapply(otu_t.m[i], as.character)
	set.seed(935234) #setting seed for same color selection
	P40 <- createPalette(16, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
	names(P40) <- sort(unique(otu_t.m$taxa))                                                               
	otu_t.m$taxa <- factor(otu_t.m$taxa, levels=unique(otu_t.m$taxa)) 
	p2 <- ggplot(data = otu_t.m, mapping = aes(x = variable, y = value, fill=taxa)) + geom_bar(stat="identity",  position="fill") + scale_fill_manual(values=P40) + theme (text = element_text(face = "bold", size = 24), axis.text.x = element_text(angle=90), legend.title=element_blank()) + ylab("sDMDMm2 Proportion") + xlab("Sample Names") + theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
	pdf(paste(prefix,"pdf",sep="."), width=12, height=10)
	print(p2)
	dev.off()
}



eucl.vec <- function(x1,x2) {
	x3 <- abs(x1 - x2)
	x4 <- x3^2
	#tp <- length(which(x2 >= 0.01))
	tp <- length(which(x2 >= 0.01 & rownames(x2) %in% rownames(x1) ==TRUE))
	recall <- tp/nrow(x1)
	precision <- tp/length(which(x2 != 0))
	list1 <- c(unlist(x3), sum(x3), sqrt(sum(x4)), tp, precision, recall)
	return(list1)
}

readMock <- function() {
	mock <- "data/Figure2/mock_distribution.txt"
	mockD <- read.table(mock, sep="\t", header=TRUE, row.names = 1)
	return(mockD)
}
