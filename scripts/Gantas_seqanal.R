countdata <- as.data.frame(seqdata) %>%
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  rename_all(str_remove, ".sorted.bam") %>%
  rename_all(str_remove, "./bam/") %>%# remove the ".bam" from the column names
  dplyr::select(sampleinfo$FileName) %>% # keep sample columns using sampleinfo%
  as.matrix()
# would be to replace 0 values by 1 read (as these ones are ignored)
#countdata[countdata==0]=1
keep <- rowSums(countdata) > 5
#keep <- rowSums(countdata) > 100
countdata <- countdata[keep,]
# for now design takes into account only
sampleinfo$dissected <- as.factor(sampleinfo$dissected)
sampleinfo$Replicate <- as.factor(sampleinfo$Replicate)
sampleinfo$dissected <- as.factor(sampleinfo$dissected)
design <- as.formula(~Replicate + Groupname )
#design <- as.formula(~dissected )
sampleinfo$Groupname <- factor(sampleinfo$Groupname,
                               levels = c("0nM", "20nM",'200nM','2000nM'))
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = sampleinfo,
                                     design = design  )
dds_lrt <- DESeq(ddsObj.raw, test="LRT", reduced = ~ Replicate)
res_LRT <- results(dds_lrt)
#write.csv(res_LRT,file='2305_revisions/to_publish/deseq2.csv')
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  filter(padj < 0.001)
data_time_deseq <-res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
data_time_deseq$symbol <- mapIds(org.Dm.eg.db,
                                 keys=data_time_deseq$gene,
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first")
#mapping ids
sig_res_LRT$symbol <- mapIds(org.Dm.eg.db,
                             keys=sig_res_LRT$gene,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
ddsObj.norm <- estimateSizeFactors(ddsObj.raw)
normalized_counts <- counts(ddsObj.norm, normalized=TRUE)
df<-as.data.frame(normalized_counts)
df <- tibble::rownames_to_column(df)
#mapping ids
df$symbol <- mapIds(org.Dm.eg.db,
                    keys=df$rowname,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
df2 <- melt(df[,-c(1)],id='symbol' )
names(df2)[c(2)]<-'FileName'
df2 = merge(df2,sampleinfo[c('FileName','Groupname','Sample')], by='FileName')
df3 = dcast(df2[,c('symbol','Groupname','value')], symbol ~ Groupname ,fun.agg = function(x) mean(x))
names(df3)[2:5]<-paste('conc',names(df3)[2:5],sep='_')
