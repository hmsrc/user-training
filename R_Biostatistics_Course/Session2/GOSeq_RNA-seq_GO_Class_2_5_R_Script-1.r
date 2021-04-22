#Install BioConductor packages goseq, GO.db, org.Hs.eg.db using the biocLite.R installation script.
# updated 31st Mar,2020
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
# Load the packages - - updated 31st Mar,2020
BiocManager::install("goseq")
BiocManager::install("GO.db")
BiocManager::install("org.Hs.eg.db")

#Load goseq library
library(goseq)

#Load GO.db library
library(GO.db)

#Load org.Hs.eg.db library
library(org.Hs.eg.db)

#Using R on the laptop              
#setwd("C:/Users/Tom/OneDrive/Documents/My Received Files/RCCG_HMS/RC_Biostatistics_Courses/Biostats_RNA-Seq/Supervised/Class_5")
setwd("C:/Users/Owner/Documents/WuXi_NextCODE_Genomics/HMS_MIT_Biostats_Course/Supervised/Class_2_5")
dir() 

#Get transcript lengths for entire background  - start by transposing the data frame before converting to a vector
#read table and name genesinput
genesinput=read.table("background_geneSymbol.txt",as.is=TRUE,header=FALSE) 

#Check 
head(genesinput) #genes symbols listed in column 1

#Transpose data frame and convert to vector
genes <- as.vector(t(genesinput)) #Please note that, if you want to have by columns you should use as.vector(as.matrix(test))-(t implicitly converts to a matrix). c(t(test)) works too.

#check
is.vector(genes) #worked

#calculate median transcript lengths
transcriptlengths = getlength(genes, "hg19", "geneSymbol")

#check
head(transcriptlengths) #worked!
tail(transcriptlengths) #worked!

#create a data frame with both vectors and with no column headings
df <- data.frame(genes,transcriptlengths)
head(df)

#write table - no column names
write.table(df, file="background_geneSymbol_transcript_length.txt", sep="\t", quote=F, col.names=F, row.names=F)

##Set up for GOseq analysis
#Run each cell list with different input file
#Name GOSeqInputFile 'GOSeqInputMDA2'
GOSeqInputMDA2 = read.table("outedgeR05_geneSymbol.txt")
#GOSeqInputMDA2 = read.table("outedgeR001_geneSymbol.txt")
#GOSeqInputMDA2 = read.table("outedgeR1e_06_geneSymbol.txt")

#Convert to a single named vector containing two pieces of gene information: Gene Symbol & 0 (normal) & 1 (abnormal expression) calls
MDA2.gene.vector=GOSeqInputMDA2[,2]
names(MDA2.gene.vector)=GOSeqInputMDA2[,1]

#Name TranscriptLengthInputFile 'GOSeqGeneSymbolTranscriptLength'
#GOSeqGeneSymbolTranscriptLength = read.table("background_geneSymbol_transcript_length.txt")
GOSeqGeneSymbolTranscriptLength = df

#Convert to a single named vector containing two pieces of gene information: Gene Symbol & Associated Transcript Lengths
transcript.lengths=GOSeqGeneSymbolTranscriptLength[,2]
names(transcript.lengths)=GOSeqGeneSymbolTranscriptLength[,1]

#Grab transcript length information and Fit the probability Weighting Function and plot the resulting fit
MDA2.pwf=nullp(MDA2.gene.vector,"hg19","geneSymbol", bias.data=transcript.lengths)

#Wallenius approximation
#The use_genes_without_cat=TRUE argument allows unannotated genes in the analysis 
MDA2.GO.wall=goseq(MDA2.pwf,"hg19","geneSymbol", use_genes_without_cat=FALSE)
#MDA2.GO.wall=goseq(MDA2.pwf,"hg19","geneSymbol", use_genes_without_cat=TRUE)

#hypergeometric 
#The use_genes_without_cat=TRUE argument allows unannotated genes in the analysis 
#MDA2.GO.wall=goseq(MDA2.pwf,"hg19","geneSymbol", use_genes_without_cat=FALSE, method="Hypergeometric")
#MDA2.GO.wall=goseq(MDA2.pwf,"hg19","geneSymbol", use_genes_without_cat=TRUE, method="Hypergeometric")

#Find GO to gene symbol matches, only need to do this once for the data set.
pwf=MDA2.gene.vector
gene2cat = getgo(names(pwf), "hg19","geneSymbol")
names(gene2cat) = names(pwf)
tmp=unlist(gene2cat,use.names=FALSE)
names(tmp)=rep(names(gene2cat),times=as.numeric(summary(gene2cat)[,1]))
cat2gene=split(names(tmp),as.vector(tmp))   

#number of genes for each GO
lens=lapply(cat2gene, length)

#add adjusted p_value
data=cbind(MDA2.GO.wall[,c(1,2)], p.adjust(MDA2.GO.wall$over_represented_pvalue,method="BH"))

#query by adjusted p-value, get a smaller list
table=data[which(data[,3]<=0.05),]

#get all terms
terms <- stack(lapply(mget(table$category, GOTERM, ifnotfound=NA), Term))

#add term
table$GO_Term <- with(table, terms$values[match(terms$ind, table$category)])

#get total number of genes for each GO
table$Gene_total <- unlist(lens[table$category])

#Get all abnormal expressed genes
abexpress=names(pwf[pwf==1]) 

#Get gene list of each GO, require the gene exist in the input file (the file with 1 and 0)
genes=lapply(table$category, function(x) intersect(unlist(cat2gene[x]), abexpress))

#Add the number of observed gene to table
table$Gene_GO_Term=unlist(lapply(genes, length))

#convert the list of genes to list of concatenated gene in string format
glist=lapply(genes, function(x) paste(x, collapse = ', '))

#All the gene string to table
table$Gene_ids=unlist(glist)

#Change Column Names
colnames(table)=c("GO_Term_Id", "pvalue", "BH_pvalue", "GO_Term", "GO_Term_Gene_Total","Gene_Count_GO_Term", "Gene_Ids")

#Reorder Column Names
table=table[,c(1,4,2,3,5,6,7)]

#Write 05 tables
write.table(table, file="edgeR_05_GOseq_GO_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_GO_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_GO_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_GO_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)

#Write 001 tables
#write.table(table, file="edgeR_001_GOseq_GO_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_GO_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_GO_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_GO_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)

#Write 1e_06 tables
#write.table(table, file="edgeR_1e_06_GOseq_GO_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_GO_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_GO_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_GO_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
