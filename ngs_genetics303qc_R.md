# HMS Research Computing: Intro to Next-Gen Sequencing Technologies Part III

This is the second practicum portion of "Intro to Next-Gen Sequencing Technologies," Spring 2018 for HMS Genetics 303qc. It builds on Part II to perform differential expression analysis of a toy mouse dataset using edgeR.

You will need to have R version > 3.0 installed.

```r
#set your working directory to where the ReadPerGene.out.tab files are located:

setwd("path/to/your/counts/files")

#install edgeR and biomaRt (only do 1x)
source('http://bioconductor.org/biocLite.R')
biocLite("edgeR")
biocLite("biomaRt")

#initialize edgeR library (do every time)
library(edgeR)

#read in files, take columns 1 (annotation) and 2 (unstranded library prep), get rid of extra alignment report lines
files=dir(pattern="*\\.tab$")
dge<-readDGE(files, columns=c(1,2), header=FALSE)
realgene <- grep("^ENS",rownames(dge))
dge <- dge[realgene,]

#set annotation
Group<-factor(c(rep("EpcamPlus", 2), rep("EpCamMinus", 2)))

colnames(dge)=c("EpcamPlusR1", "EpcamPlusR2", "EpcamMinusR1", "EpcamMinusR2")

#make the edgeR object, combining counts and group
y<-DGEList(counts=dge, group=Group)

#drop annotation that does not have 3 cpm in at least 2 samples. ~ (50%) filtering
keep <- rowSums(cpm(y)>3) >= 2
filter <- y[keep,,keep.lib.sizes=FALSE]
dim(filter) #12602

#TMM normalization
norm<-calcNormFactors(filter)
norm$samples

#multidimensional scaling plot
#replicates should separate along X axis first, then y axis
png("MDS_run1.png")
plotMDS(norm, col=as.numeric(Group), main="MDS")
dev.off()

#to display as dots
png("MDS_dots_run1.png")
plotMDS(norm, col=as.numeric(Group), main="MDS", pch=18)
dev.off()

#estimate dispersion, plot overall BCV averaged over all genes
common <- estimateCommonDisp(norm, verbose=TRUE)
#Disp = 0.10518 , BCV = 0.3243 

#estimate gene-specific dispersion, plot
tag <- estimateTagwiseDisp(common)

png("BCV_run1.png")
plotBCV(tag, main="Biological Coeffient of Variation")
dev.off()

#fit the negative binomial, with Plus over Minus (2 over 1)
et<-exactTest(tag, pair=2:1)

#perform the decideTest, those making q<0.05
summary(de <- decideTestsDGE(et))

#pull only significant genes
detags <- rownames(et)[as.logical(de)]

#create significant gnenes subset normalized counts 
mydf<-norm[which(row.names(norm) %in% detags),]

#tell edgeR how many genes to pull metrics from, based on length of significant genes
lendetags<-length(detags)
tab<-topTags(et, n=lendetags)

#create data frame pulling cpm + metrics for only BH results
annot<-merge(cpm(mydf, normalized.lib.sizes=TRUE), tab, by="row.names")

#write out final norm cpm + stats table
filename <- paste("exact.test", ".run1.txt", sep="")
write.table(annot, filename, sep="\t")

#smear MDA-like plot

pngname <- paste("exact.test", ".run1.plotSmear.png")
png(pngname)
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()
```

### biomaRt Annotation pull

Creates an additional output file with gene symbol information added (takes more time, and sometimes biomaRt is down)

```r
library(biomaRt)

#set mart to mouse
mart = useEnsembl("ENSEMBL_MART_ENSEMBL")
mart=useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")

#read in file ending in .run1
files<-dir(pattern="*\\.run1.txt$")

input<-files[1] #file 1
df<-read.table(input, header=T, sep="\t") #read table
ensembl_ids<-df[,1] #get ensembl ids, first column
ids<-NULL
#get gene name and ensembl gene for each ensembl id
ids<-getBM(attributes=c("external_gene_name", "ensembl_gene_id"), filters="ensembl_gene_id", values=ensembl_ids, mart=mart) 
#mrege together input with gene symbols, drops unannotated samples
output<-merge(ids, df, by.x="ensembl_gene_id", by.y="Row.names")
colnames(output)[1:2]<-c("Ensembl", "GeneID")
filename<-gsub(".txt", ".ids.txt", input)
write.table(output, file=filename, sep="\t", quote=FALSE, row.names=FALSE)
