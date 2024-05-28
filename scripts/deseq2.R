library(data.table)
library(dplyr)
library(tibble)
library(ggplot2)
counts <- fread(file ="rsem.merged.gene_counts.tsv", header=TRUE)
gcounts <- counts[,-2]

dt <- as.data.frame(gcounts) %>%
  column_to_rownames("gene_id")  # turn the geneid column into rownames

dt[,1:9] <- lapply(dt[,1:9], as.integer)

mcountdata <- as.matrix(dt)  
  
sample_info <- read.delim(file = "RN23010-sampleinfo_full.txt")
#new_info <- sample_info %>%  filter(!row_number() %in% c(1, 2, 3))
info<-sample_info %>% column_to_rownames("Sample")
has_rownames(dt)

#I don't understand how this command works
all(rownames(info) == colnames(mcountdata))
mcountdata_s <- mcountdata[, rownames(info)]
all(rownames(info) == colnames(mcountdata_s))

#why do we need to do this?
info$Group <- as.factor(info$Group)
info$Repeat <- as.factor(info$Repeat)
#Do I need to do that
info$Group <- factor(info$Group,levels = c("control","sal-brk", "ubi-brk"))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = mcountdata_s,
                              colData = info,
#The object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis is the DESeqDataSet, which will usually be represented in the code here as an object dds.
design = ~Group)
dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$log2FoldChange),]
summary(res)
plotMA(res, ylim=c(-4,4))
plotCounts(dds, gene="FBgn0024250", intgroup="Group") #Brk counts across samples
plotCounts(dds, gene="FBgn0000179", intgroup="Group") #omb norm counts
plotCounts(dds, gene="FBgn0003525", intgroup="Group") #Sal
plotCounts(dds, gene="FBgn0025360", intgroup="Group") #optix
plotCounts(dds, gene="FBgn0011706", intgroup="Group") #reaper
plotCounts(dds, gene="FBgn0003525", intgroup="Group") #string
plotCounts(dds, gene="FBgn0262656", intgroup="Group") #myc
plotCounts(dds, gene="FBgn0010770", intgroup="Group") #ppan‹‹

rld = rlog(dds)
z <- plotPCA(rld, intgroup = "Group")
z + geom_text_repel(aes(x=PC1, y=PC2, label=rownames(info)))+
  theme_bw()
