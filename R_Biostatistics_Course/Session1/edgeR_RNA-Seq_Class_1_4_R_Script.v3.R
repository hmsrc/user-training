##edgeR Differential Expression analyses of RNA-seq data 
##Based edgeR R Tutorial by Sean Ruddy, 2011
##Based on Edger:differential expression analysis of digital gene expression data User's Guide. Yunshun Chen, Davis McCarthy, Mark Robinson, Gordon K. Smyth, 2014.
##Citation: Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26, pp. -1.

##Working on O2
#Login to 02, with 2 factor outside of HMS private
ssh -XY ecommonas@o2.hms.harvard.edu
#Work on )2 Interactively
#Interactive Jobs - Max of 12 hours
srun --pty -p interactive -t 0-12:00 --mem 8G bash
#You will then be assigned a Compute Node, which will allow you to work interactively for up to 12 hours.

#Interactive Jobs - Max of 12 hours & 32G of memory
srun --pty -p interactive -t 0-12:00 --mem 32G bash

#Setting up a personal R library on O2 
mfk8@login00:~$ mkdir -p ~/R/library
mfk8@login00:~$ echo 'R_LIBS_USER="~/R/library"' >  $HOME/.Renviron
mrk8@login00:~$ export R_LIBS_USER="/home/user123/R/library"

#set working directory - Example
cd /your/working/Directory/

#Directly Set R Version
module load gcc/6.2.0 R/3.5.1-extra

#Launch R
R

##Or working on a Mac
##setwd("/your/working/Directory/")
dir()

####Class 1-4 - Assessment of Differential Expression of RNA-seq Data with edgeR####
#Install edgeR - updated 31st Mar,2020
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# updated 31st Mar,2020
BiocManager::install("edgeR")

#Load the edgeR package
library(edgeR)

#set working directory
#setwd("C:/Users/Owner/Documents/RCCG_HMS/RC_Biostatistics_Courses/Biostats_RNA-Seq/Supervised/Class_4")
#setwd("C:/Users/Tom/OneDrive/Documents/My Received Files/RCCG_HMS/RC_Biostatistics_Courses/Biostats_RNA-Seq/Supervised/Class_4")
setwd("C:/Users/Owner/Documents/WuXi_NextCODE_Genomics/HMS_MIT_Biostats_Course/Supervised/Class_1_4")
dir()

#read in entire counts data .txt file
x <- read.delim("TNBC10vNormal10_Raw_2.txt", row.names=1, stringsAsFactors=FALSE)
head(x)
tail(x)
dim(x)

#create groups
group <- c(rep("TNBC", 10), rep("Normal", 10))

#Check
group

#Building the edgeR Object
#create DGEList object to work on, where are the columns are the counts
#Put the counts and other information into a DGEList object
#colnames(y) <- targets$Label
y <- DGEList(counts=x[,1:20], group=group)
dim(y)
#[1] 20531    20


#Filter out low count reads since it would be impossible to detect differential expression. The method used in the
#edgeR vignette is to keep only those genes that have at least 1 read per million in at least 10 samples. Therefore, Filter out very lowly expressed tags, 
#keeping genes that are expressed at a reasonable level in at least one condition. Since the group sizes are 10, keep genes that achieve at least 
#one count per million (cpm) in 10 samples
#Once this filtering is done, calculate the normalization factors which correct for the different compositions of the samples. The effective
#library sizes are then the product of the actual library sizes and these factors.
keep <- rowSums(cpm(y)>1) >= 10
y <- y[keep,]
dim(y)
#[1] 14332    20

#Re-compute the library sizes
y$samples$lib.size <- colSums(y$counts)
#Check
y$samples

##TMM Normalization
#Compute effective library sizes using Trimmed Mean of M Values (TMM) normalization
#This accounts for differences in transcriptome sizes
y <- calcNormFactors(y)
y$samples

#effective library sizes
y$samples$lib.size*y$samples$norm.factors

##Data exploration
#Generate a Multi-dimensional Scaling (MDS) plot
#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions.
#An MDS plot show distances, in terms of biological coefficient of variation (BCV), between samples
#Dimension 1 separates the Normal from TNBC samples. The MDS plot shows the replicates are consistent. 
#This is a good indication edgeR analysis will produce adequate DE genes.
#To view MDS plot
plotMDS(y, col=as.numeric(as.factor(group)), main ="MDS Plot for TMM Normalized Count Data")
?cmdscale

#Write plot as a pdf
pdf("MDS_Plot.pdf" , width = 11 , height = 7) # in inches
plotMDS(y, col=as.numeric(as.factor(group)), main ="MDS Plot for TMM Normalized Count Data")
dev.off() # this tells [R] to close and stop writing to the pdf

##Estimating common and tagwise dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes with the qCML method
y <- estimateCommonDisp(y, verbose=TRUE)
#Disp = 0.39294 , BCV = 0.6268
#BCV (square root of the common dispersion) is *63%. BCV for a laboratory experiment with cell lines or model organisms is typically 10-15%.

#Now estimate gene-specific (tagwise) dispersions with the qCML method
y <- estimateTagwiseDisp(y)
summary(y$tagwise.dispersion)
#Min.    1st Qu. Median  Mean    3rd Qu. Max. 
#0.1097  0.1948  0.2763  0.3864  0.4544  4.4520 

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y, main = "Plot of Estimated Dispersions")

#Write plot as a pdf
pdf("Est_Dispersions_Plot.pdf" , width = 11 , height = 7) # in inches
plotBCV(y, main = "Plot of Estimated Dispersions")
dev.off() # this tells [R] to close and stop writing to the pdf

##Differential expression
#Compute NB exact genewise tests for differential expression between TNBC and Normal samples
et <- exactTest(y)
options(digits = 3) #print only 3 digits
top <- topTags(et)
top

#Check the individual cpm values for the top genes
#Check up/down expression to TNBC/Normal *VERY IMPORTANT*
cpm(y)[rownames(top), ]

#print the total number of DE genes at 5% FDR
summary(de <- decideTestsDGE(et))
#   [,1]
#-1 2290
#0  9019
#1  3023 
##14,332 total genes after filtering for low expression
##Of the 5313 genes identified as DE, 3023 are up-regulated in TNBC samples and 2290 are down-regulated.

#Plot the log-fold-changes with a smear plot, highlighting the DE genes
#Blue lines indicate 2-fold changes
#to view Smear Plot
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags, main = "Smear Plot")
abline(h=c(-1, 1), col="blue")

#Write plot as a pdf
pdf("Smear_Plot.pdf" , width = 7 , height = 7) # in inches
plotSmear(et, de.tags=detags, main = "Smear Plot")
abline(h=c(-1, 1), col="blue")
dev.off() # this tells [R] to close and stop writing to the pdf

#Take a look at the smear plot function
?plotSmear

#Create an object called "resultsTbl.tagwise" containing the four output parameters: logFC, logCPM, PValue,FDR
resultsTbl.tagwise <- topTags(et, n = nrow(et$table))$table

#Check resultsTbl.tagwise
head(resultsTbl.tagwise)
dim(resultsTbl.tagwise)

#Generate objects according to three significance levels for tagwise dispersion estimates: 0.05, 0.001, 1e-06
#These will be used for GOseq analyses in classes 5 and 6.
#Names/IDs of DE genes
de.genes.tagwise_0.05 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 0.05]
de.genes.tagwise_0.001 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 0.001]
de.genes.tagwise_1e_06 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 1e-06]

#Check number of genes at each significance level
length(de.genes.tagwise_0.05) #5313 DE genes
length(de.genes.tagwise_0.001) #2398 DE genes
length(de.genes.tagwise_1e_06) #944 DE genes

#Check percentage of total genes
length(de.genes.tagwise_0.05)/nrow(resultsTbl.tagwise)*100 #37.07089
length(de.genes.tagwise_0.001)/nrow(resultsTbl.tagwise)*100 #16.73179
length(de.genes.tagwise_1e_06)/nrow(resultsTbl.tagwise)*100 #6.586659

#Check Up/Down regulated summaries for tagwise results
summary(de <- decideTestsDGE(et, p.value=0.001))
#    [,1] 
#-1  1179
#0  11934
#1   1219

summary(de <- decideTestsDGE(et, p.value=1e-06))
#    [,1] 
#-1   535
#0  13388
#1    409

#Check de.genes.tagwise_1e_06
head(de.genes.tagwise_1e_06) #Note de.genes.tagwise_1e_06 is a vector of 1e06 genes

#Double check
is.vector(de.genes.tagwise_1e_06)
#[1] TRUE

##Format Results and count data form DGEList object(y)
colnames(resultsTbl.tagwise) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise <- match(rownames(resultsTbl.tagwise), rownames(y$counts))

#Check for all 14332 genes
length(wh.rows.tagwise)

#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
#combResults.tagwise <- cbind(resultsTbl.tagwise,
#                             "Tgw.Disp" = y$tagwise.dispersion[wh.rows.tagwise],
#                             "UpDown.Tgw" = decideTestsDGE(et, p.value = 0.05)[wh.rows.tagwise],
#                             y$counts[wh.rows.tagwise,])
                             
#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
combResults.tagwise <- cbind(resultsTbl.tagwise,                             "Tgw.Disp" = y$tagwise.dispersion[wh.rows.tagwise],
                            "UpDown.Tgw" = edgeR::decideTestsDGE(et, p.value = 0.05)@.Data[,1][wh.rows.tagwise],                                y$counts[wh.rows.tagwise,])
                               

#Check combResults.tagwise 
head(combResults.tagwise)
dim(combResults.tagwise)

#Create significance level objects from original counts data (x) and FDR vectors
de.tag_05<-x[which(row.names(x)%in%(de.genes.tagwise_0.05)),]
de.tag_001<-x[which(row.names(x)%in%(de.genes.tagwise_0.001)),]
de.tag_1e_06<-x[which(row.names(x)%in%(de.genes.tagwise_1e_06)),]

#Check
head(de.tag_05)
dim(de.tag_05)
head(de.tag_001)
dim(de.tag_001)
head(de.tag_1e_06)
dim(de.tag_1e_06)

#Check again
head(combResults.tagwise)
head(resultsTbl.tagwise)

##Merge all output into FDR specific objects
combResults.tagwise_05<-merge(combResults.tagwise [,1:6], de.tag_05, by="row.names")
head(combResults.tagwise_05) #Almost there. Need to rename first column and Ensembl and Entrez columns to front. Ensembl and Entrez columns are there as a double check

#Rename column 1 as "HGNC"
colnames(combResults.tagwise_05)[1] ="HGNC"
head(combResults.tagwise_05)
dim(combResults.tagwise_05)

#Create edgeR_05 output object and reorder columns 
edgeR_05=cbind(combResults.tagwise_05[,28:29],combResults.tagwise_05[,1:27])
head(edgeR_05)
dim(edgeR_05)

#Repeat for FDR of 0.001
combResults.tagwise_001<-merge(combResults.tagwise [,1:6], de.tag_001, by="row.names")
head(combResults.tagwise_001)

#Rename column 1 as "HGNC"
colnames(combResults.tagwise_001)[1] ="HGNC"
head(combResults.tagwise_001)
dim(combResults.tagwise_001)

#Create edgeR_001 output object and reorder columns
edgeR_001=cbind(combResults.tagwise_001[,28:29],combResults.tagwise_001[,1:27])
head(edgeR_001)
dim(edgeR_001)

#Repeat for FDR of 1e06
combResults.tagwise_1e_06<-merge(combResults.tagwise [,1:6], de.tag_1e_06, by="row.names")
head(combResults.tagwise_1e_06)

#Rename column 1 as "HGNC"
colnames(combResults.tagwise_1e_06)[1] ="HGNC"
head(combResults.tagwise_1e_06)
dim(combResults.tagwise_1e_06)

#Create edgeR_1e_06 output object and reorder columns
edgeR_1e_06=cbind(combResults.tagwise_1e_06[,28:29],combResults.tagwise_1e_06[,1:27])
head(edgeR_1e_06)
dim(edgeR_1e_06)

#Sort results by FDR as a check for row names.
edgeR_05_out<-edgeR_05[order(edgeR_05[,7],decreasing=FALSE),]
head(edgeR_05_out)

edgeR_001_out<-edgeR_001[order(edgeR_001[,7],decreasing=FALSE),]
head(edgeR_001_out)

edgeR_1e_06_out<-edgeR_1e_06[order(edgeR_1e_06[,7],decreasing=FALSE),]
head(edgeR_1e_06_out)

#Write all three edgeR FDR specific objects to .txt files - row.names=FALSE to maintain appropriate column order 
write.table(edgeR_05_out, file="edgeR_tagwise_05_output.txt", row.names=FALSE, sep="\t")
write.table(edgeR_001_out, file="edgeR_tagwise_001_output.txt", row.names=FALSE, sep="\t")
write.table(edgeR_1e_06_out, file="edgeR_tagwise_1e_06_output.txt", row.names=FALSE, sep="\t")



