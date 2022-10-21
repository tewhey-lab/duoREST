
#R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"


library(DESeq2)
df <- read.delim('RE1KO_BRBseq_output.dge.umis.detailed.txt', sep='\t', header=TRUE)
df1 <- df[grepl('ENSG', df$Gene_id),]
df1 <- df1[,grepl('_r', colnames(df1))]
df1 <- df1[,order(names(df1))]
NTs <- df1[, grepl('NT', colnames(df1))]
coldata <- read.delim('fig_5_BRBseq_colData.txt', sep='\t',header=TRUE, row.names=1)


#### running target1 ####
VC1s <- cbind(NTs, df1[, grepl('target1', colnames(df1))])
coldataVC1 <- rbind(coldata[grepl('NT', rownames(coldata)),], coldata[grepl('target1', rownames(coldata)),])
ddsVC1 <- DESeqDataSetFromMatrix(countData = VC1s, colData=coldataVC1, design = ~ target)
ddsVC1$target <- relevel(ddsVC1$target, 'NT')
resVC1 <- DESeq(ddsVC1)
resultsVC1 <- results(resVC1)
ddsVC1 <- estimateSizeFactors(ddsVC1)
countsVC1 <- counts(ddsVC1, normalized=TRUE)

df2 <- df[grepl('ENSG', df$Gene_id),]

countsVC1 <- cbind(df2$Gene_name, countsVC1)
colnames(countsVC1)[1] <- "Gene_name"
write.table(countsVC1, 'norm_count.target1.txt', sep='\t', col.names=TRUE, row.names = FALSE)

resultsVC1mtrx <- cbind(df2$Gene_name, resultsVC1)
colnames(resultsVC1mtrx)[1] <- "Gene_name"
write.table(resultsVC1mtrx, 'DESeq.target1.txt', sep='\t', col.names=TRUE, row.names = FALSE)

#### running target2 ####

VC2s <- cbind(NTs, df1[, grepl('target2', colnames(df1))])
coldataVC2 <- rbind(coldata[grepl('NT', rownames(coldata)),], coldata[grepl('target2', rownames(coldata)),])
ddsVC2 <- DESeqDataSetFromMatrix(countData = VC2s, colData=coldataVC2, design = ~ target)
ddsVC2$target <- relevel(ddsVC2$target, 'NT')
resVC2 <- DESeq(ddsVC2)
resultsVC2 <- results(resVC2)
ddsVC2 <- estimateSizeFactors(ddsVC2)
countsVC2 <- counts(ddsVC2, normalized=TRUE)

countsVC2 <- cbind(df2$Gene_name, countsVC2)
colnames(countsVC2)[1] <- "Gene_name"
write.table(countsVC2, 'norm_count.target2.txt', sep='\t', col.names=TRUE, row.names = FALSE)

resultsVC2mtrx <- cbind(df2$Gene_name, resultsVC2)
colnames(resultsVC2mtrx)[1] <- "Gene_name"
write.table(resultsVC2mtrx, 'DESeq.target2.txt', sep='\t', col.names=TRUE, row.names = FALSE)


#### running target3 ####

VC3s <- cbind(NTs, df1[, grepl('target3', colnames(df1))])
coldataVC3 <- rbind(coldata[grepl('NT', rownames(coldata)),], coldata[grepl('target3', rownames(coldata)),])
ddsVC3 <- DESeqDataSetFromMatrix(countData = VC3s, colData=coldataVC3, design = ~ target)
ddsVC3$target <- relevel(ddsVC3$target, 'NT')
resVC3 <- DESeq(ddsVC3)
resultsVC3 <- results(resVC3)
ddsVC3 <- estimateSizeFactors(ddsVC3)
countsVC3 <- counts(ddsVC3, normalized=TRUE)

countsVC3 <- cbind(df2$Gene_name, countsVC3)
colnames(countsVC3)[1] <- "Gene_name"
write.table(countsVC3, 'norm_count.target3.txt', sep='\t', col.names=TRUE, row.names = FALSE)

resultsVC3mtrx <- cbind(df2$Gene_name, resultsVC3)
colnames(resultsVC3mtrx)[1] <- "Gene_name"
write.table(resultsVC3mtrx, 'DESeq.target3.txt', sep='\t', col.names=TRUE, row.names = FALSE)


#### running target4 ####

VC4s <- cbind(NTs, df1[, grepl('target4', colnames(df1))])
coldataVC4 <- rbind(coldata[grepl('NT', rownames(coldata)),], coldata[grepl('target4', rownames(coldata)),])
ddsVC4 <- DESeqDataSetFromMatrix(countData = VC4s, colData=coldataVC4, design = ~ target)
ddsVC4$target <- relevel(ddsVC4$target, 'NT')
resVC4 <- DESeq(ddsVC4)
resultsVC4 <- results(resVC4)
ddsVC4 <- estimateSizeFactors(ddsVC4)
countsVC4 <- counts(ddsVC4, normalized=TRUE)

countsVC4 <- cbind(df2$Gene_name, countsVC4)
colnames(countsVC4)[1] <- "Gene_name"
write.table(countsVC4, 'norm_count.target4.txt', sep='\t', col.names=TRUE, row.names = FALSE)

resultsVC4mtrx <- cbind(df2$Gene_name, resultsVC4)
colnames(resultsVC4mtrx)[1] <- "Gene_name"
write.table(resultsVC4mtrx, 'DESeq.target4.txt', sep='\t', col.names=TRUE, row.names = FALSE)


#### running target5 ####

VC5s <- cbind(NTs, df1[, grepl('target5', colnames(df1))])
coldataVC5 <- rbind(coldata[grepl('NT', rownames(coldata)),], coldata[grepl('target5', rownames(coldata)),])
ddsVC5 <- DESeqDataSetFromMatrix(countData = VC5s, colData=coldataVC5, design = ~ target)
ddsVC5$target <- relevel(ddsVC5$target, 'NT')
resVC5 <- DESeq(ddsVC5)
resultsVC5 <- results(resVC5)
ddsVC5 <- estimateSizeFactors(ddsVC5)
countsVC5 <- counts(ddsVC5, normalized=TRUE)

countsVC5 <- cbind(df2$Gene_name, countsVC5)
colnames(countsVC5)[1] <- "Gene_name"
write.table(countsVC5, 'norm_count.target5.txt', sep='\t', col.names=TRUE, row.names = FALSE)

resultsVC5mtrx <- cbind(df2$Gene_name, resultsVC5)
colnames(resultsVC5mtrx)[1] <- "Gene_name"
write.table(resultsVC5mtrx, 'DESeq.target5.txt', sep='\t', col.names=TRUE, row.names = FALSE)


