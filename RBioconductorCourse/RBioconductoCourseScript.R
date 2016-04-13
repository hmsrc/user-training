####
#R code for R/Bioconductor Course - HMS Research Computing - 04/14/2016
#kristina_holton@hms.harvard.edu
####

#setwd to where you saved Rcoursetestdata1.csv and Rcoursetestdata2.csv
setwd("location/of/data")

#read in Rcoursetestdata1.csv as a data frame with header and row names
mydf <- read.table("Rcoursetestdata1.csv", header=TRUE, row.names=1, sep=",")

#print mydf
mydf

#get summary statistics on mydf
summary(mydf)

#turn mydf into a new matrix
mymatrix <- as.matrix(mydf)

#transpose mymatrix - only matrices can be transposed
myTmatrix<- t(mymatrix) 

#turn myTmatrix back into a data frame
myTdf <- as.data.frame(myTmatrix) 

#print myTdf
myTdf

#get summary statistics on myTdf
summary(myTdf)

####
#Statistics

#perform a t test over TNBC (columns 1-3) versus Normal (columns 4-6).  
t.test(mydf[,1:3], mydf[,4:6])

#perform a Wilcoxon test over TNBC (columns 1-3) versus Normal (columns 4-6), requires numerical matrix
wilcox.test(mymatrix[,1:3], mymatrix[,4:6])

#dumby test of a linear model using TNBC1 as dependent, Normal1 as independent variables
mymodel  <- lm(TNBC1 ~ Normal1, data=mydf)
anova(mymodel)

####
#Plotting

#create a boxplot over samples from mydf
boxplot(mydf, main="Boxplot of My Samples", xlab="Sample", ylab="Gene Values", col=c("red", "orange", "yellow", "green", "blue", "purple"))

#create a boxplot over genes fromm myTdf
boxplot(myTdf, main="Boxplot of My Genes", xlab="Genes", ylab="Gene Values")

#create a barplot of genes from matrix
barplot(mymatrix, main="Barplot of My Genes", xlab="Gene", ylab="Gene Values", col=c("blue"))

#create a histogram of mymatrix values
hist(mymatrix, main="Histogram of My Matrix", col=c("red", "orange", "yellow", "green", "blue", "purple"))

#for hierarchical clustering: create distance matrix
d <- dist(myTmatrix) #takes matrix as input
hc <- hclust(d) #creates hcl matrix
plot(hc)

#create a line plot
#start a line plot with ENSG00000008988 "b" means both points and lines
plot(myTdf$ENSG00000008988, type="b", col="green",  ylim=c(10000,150000), main="Gene Values Over Samples", xlab="Sample", ylab="Gene Values")
#add a new line for ENSG00000009307 "lines" adds a line to the current plot
lines(myTdf$ENSG00000009307, type="b", col = "blue")
#add a new line for ENSG00000019582
lines(myTdf$ENSG00000019582, type="b", col="red")
#create a legend
legend(5, 140000, #positions x, y
c("TPS20", "CSDE1", "CD74"), #line names
lty=c(1,1), #specifies lines
lwd=c(2.5,2.5), #specifies line width
col=c("green", "blue", "red") #add colors
) #ends legend


####
#Count Data

#read in Rcoursetestdata2.csv counts table as a data frame
mytable<-read.table("Rcoursetestdata2.csv", header=T, row.names=1, sep=",")

#print mytable
mytable

#summary statistic of mytable
summary(mytable)

#create a table of gender vs smoke
genderVsmoke <- table(mytable$gender, mytable$smoke)

#print genderVsmoke
genderVsmoke

#do a summary chi square test on genderVsmoke
summary(genderVsmoke)

#do a mosaic plot of genderVsmoke
mosaicplot(genderVsmoke)

####
#Bonus Heatmap
#uses "heatmap.2" from package "gplots" and color ramped palette from "RColorBrewer"

#install packages and call libraries

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

#define a Brewer palette
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)

#define color breaks based on data
colors = c(seq(0,10000,length=100),seq(10001,50000,length=100),seq(50001,150000,length=100))

#set input data matrix
mydata <- mymatrix

heatmap.2(mydata,
main="TNBC + Normals Heatmap", #Title
xlab= "Samples", # x axis label
ylab= "Genes", # y axis label
density.info="none", #turns off density plot inside color legend
trace="none", #turns off trace lines inside heatmap
col=my_palette, #uses color palette defined above
breaks=colors, #uses color breaks defined above
hclustfun=function(x) hclust(x,method="complete"), #specifies type of linkage
distfun=function(x) dist(x,method="euclidean"), #specifies distance metric
dendrogram="column", #draws a column dendrogram 
srtCol=0, #horizontal column names (0 degrees)
offsetCol=2,#column name spacing
adjCol=c(0.5,1), #spacing between column names
labRow=FALSE #turn off gene labels (too many)
) #close plot
