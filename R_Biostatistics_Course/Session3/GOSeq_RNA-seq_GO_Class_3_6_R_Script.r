#Install BioConductor packages goseq, GO.db, org.Hs.eg.db using the biocLite.R installation script.
source("http://bioconductor.org/biocLite.R")


#biocLite("gage")
biocLite("pathview")

#Load goseq library
library(goseq)

#Load GO.db library
library(KEGG.db)

#Load org.Hs.eg.db library
library(org.Hs.eg.db)

#Load gage library
#library(gage)

#Load pathview library
library(pathview)

#Using R on the laptop              
#setwd("C:/Users/Owner/Documents/RCCG_HMS/RC_Biostatistics_Courses/Biostats_RNA-Seq/Supervised/Class_6")
#setwd("C:/Users/Tom/OneDrive/Documents/My Received Files/RCCG_HMS/RC_Biostatistics_Courses/Biostats_RNA-Seq/Supervised/Class_6")
setwd("C:/Users/Owner/Documents/WuXi_NextCODE_Genomics/HMS_MIT_Biostats_Course/Supervised/Class_3_6")
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
MDA2.KEGG.wall=goseq(MDA2.pwf,"hg19","geneSymbol", test.cats="KEGG", use_genes_without_cat=FALSE)
#MDA2.KEGG.wall=goseq(MDA2.pwf,"hg19","geneSymbol", test.cats="KEGG", use_genes_without_cat=TRUE)

#hypergeometric 
#The use_genes_without_cat=TRUE argument allows unannotated genes in the analysis 
#MDA2.KEGG.wall=goseq(MDA2.pwf,"hg19","geneSymbol", test.cats="KEGG", use_genes_without_cat=FALSE, method="Hypergeometric")
#MDA2.KEGG.wall=goseq(MDA2.pwf,"hg19","geneSymbol", test.cats="KEGG", use_genes_without_cat=TRUE, method="Hypergeometric")

#Find KEGG to gene symbol matches, only need to do this once for the data set.
pwf=MDA2.gene.vector
gene2cat = getgo(names(pwf), "hg19","geneSymbol", fetch.cats=c("KEGG"))
names(gene2cat) = names(pwf)
tmp=unlist(gene2cat,use.names=FALSE)
names(tmp)=rep(names(gene2cat),times=as.numeric(summary(gene2cat)[,1]))
cat2gene=split(names(tmp),as.vector(tmp))   

#number of genes for each KEGG
lens=lapply(cat2gene, length)

#add adjusted p_value
data=cbind(MDA2.KEGG.wall[,c(1,2)], p.adjust(MDA2.KEGG.wall$over_represented_pvalue,method="BH"))

#query by adjusted p-value, get a smaller list
table=data[which(data[,3]<=0.05),]

#get and add all terms
table$KEGG_PATH_NAME <- unlist(mget(table$category, KEGGPATHID2NAME, ifnotfound=NA))

#get total number of genes for each KEGG
table$Gene_total <- unlist(lens[table$category])

#Get all abnormal expressed genes
abexpress=names(pwf[pwf==1]) 

#Get gene list of each KEGG, require the gene exist in the input file (the file with 1 and 0)
genes=lapply(table$category, function(x) intersect(unlist(cat2gene[x]), abexpress))

#Add the number of observed gene to table
table$SNV_Gene_KEGG_Term=unlist(lapply(genes, length))

#Convert the list of genes to list of concatenated gene in string format
glist=lapply(genes, function(x) paste(x, collapse = ', '))

#Add all the gene string to table
table$Gene_ids=unlist(glist)

#Change Column Names
colnames(table)=c("KEGG_Id", "pvalue", "BH_pvalue", "KEGG_Term", "KEGG_Term_Gene_Total","Gene_Count_KEGG_Term", "Gene_Ids")

#Reorder Column Names
table=table[,c(1,4,2,3,5,6,7)]

#Write 05 tables
write.table(table, file="edgeR_05_GOseq_KEGG_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_KEGG_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_KEGG_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_05_GOseq_KEGG_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)

#Write 001 tables
#write.table(table, file="edgeR_001_GOseq_KEGG_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_KEGG_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_KEGG_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_001_GOseq_KEGG_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)

#Write 1e_06 tables
#write.table(table, file="edgeR_1e_06_GOseq_KEGG_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_KEGG_hypergeo_restrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_KEGG_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)
#write.table(table, file="edgeR_1e_06_GOseq_KEGG_hypergeo_unrestrict.txt", sep="\t", quote=F, col.names=T, row.names=F)

##Print enriched KEGG pathway with Pathview
#Read in the Entrez Gene Ids from the 0.05, 0.001, and 1e-06 analyses as data frames (read.table)
#Grab enriched KEGG pathway Ids from column 1 of the table "table" above
#The apply loop iterates over all KEGG pathways in the table
entrezids<-read.table("edgeR_tagwise_05_out_Entrez.txt")
pv.out.list<-sapply(table[,1], function(pid) pathview(gene.data=entrezids, pathway.id= pid, species = "hsa"))

