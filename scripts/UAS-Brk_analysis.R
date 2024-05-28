library(DESeq2)
library(tidyverse)
library(dplyr)

#uploading data
ubi_brk <- read_tsv("MIK5655A5.featureCounts", comment="#")
sal_brk <- read_tsv("MIK5655A4.featureCounts", comment="#")
control <- read_tsv("MIK5655A3.featureCounts", comment="#")

#changing a bulky column name
colnames(control)[7] <- "control"
colnames(ubi_brk)[7] <- "ubi_brk"
colnames(sal_brk)[7] <- "sal_brk"

#transforming table into dataframe format
control_d <-as.data.frame(control)
ubi_brk_d <-as.data.frame(ubi_brk)
sal_brk_d <-as.data.frame(sal_brk)

#merge the data in one table
mer<-merge(control_d, ubi_brk_d)
mcounts <-merge(mer, sal_brk_d)

countdata<- mcounts %>% 
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  dplyr::select(6:8) %>% # keep sample columns; rownames doesnt exist
#as a column
  as.matrix()
#saving the countdata
write.table(countdata, "countdata.tsv", sep = "\t")
#write.table(countdata, "countdata.txt", col.names=NA)

#We will keep all genes where the total number of reads across all samples 
#is greater than 5. from 23k rows to 13k
dim(countdata)
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]
dim(countdata)

#comparing library sizes across samples
librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=60e6, lty=2)

# Get log2 counts
logcounts <- log2(countdata + 1)

# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2)
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(as.matrix(logcounts)), col="blue")


#Use the DESeq2 function rlog to transform the count data. 
#This function also normalises for library size.
rlogcounts <- rlog(countdata)
write.csv(rlogcounts, "Rlogcounts.csv")
# Check distributions of samples using boxplots
boxplot(rlogcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")


#Performing PCA analysis on top of rlog we got
library(ggfortify) #this lib can recognise outputs of stat analysis and plot 
#them automatically
library(ggrepel)
# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
 autoplot(pcDat, size=4)+
  geom_text_repel(aes(x=PC1, y=PC2, label=colnames(countdata)), box.padding = 0.8)+
   theme_bw()

#Convert counts to Deseq objects
 # create the design formula
 design <- as.formula(~ colnames.countdata.)
 # create the DESeqDataSet object
 ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = DataFrame(colnames(countdata)),
                                  design = design)
 ## converting counts to integer mode
 ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
 ## design formula are characters, converting to factors
 
 
 # Apply normalisation to DDS object
 ddsObj <- estimateSizeFactors(ddsObj)
#Looking at resulting normalization factors
 sizeFactors(ddsObj)
 
 save(countdata, file="preprocessing.RData")
 
 
 #Creating the design model formula
 design <- as.formula(~ colnames.countdata.)
 modelMatrix <- model.matrix(design, data = colData)
 modelMatrix
 
 
 #Let’s plot a PCA from vst transformed data. Can you anticipate if the interaction term will be important?
 
 vstcounts <- vst(ddsObj, blind=TRUE)
 plotPCA(vstcounts, intgroup=c("colnames.countdata."))+
   geom_text_repel(aes(x=PC1, y=PC2, label=colnames(countdata)), box.padding = 0.8)+
   theme_bw()
 #################################################################################
 #… then estimate dispersion …
 
 ddsObj <- estimateDispersions(ddsObj)
#The design matrix has the same number of samples and coefficients to fit,
 # so estimation of dispersion is not possible. Treating samples
 # as replicates was deprecated in v1.20 and no longer supported since v1.22. 
 # 
 
 #… finally, apply Negative Binomial GLM fitting and calculate Wald statistics
 ddsObj <- nbinomWaldTest(ddsObj)
 # Error in nbinomWaldTest(ddsObj) : 
 #   testing requires dispersion estimates, first call estimateDispersions()
 
 # Run DESeq
 ddsObj <- DESeq(ddsObj)
 #################################################################################
 
 
 #Generate a results table
 res <- results(ddsObj, alpha=0.05)
 can I run deseq if i dint have replicates
 head(res)
 
 
 
 