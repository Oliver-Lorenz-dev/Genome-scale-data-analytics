library('DESeq2')
#Task 7.1
load(file = "/home/student/practical_2/sim_data/Ctrachomatis_counts.RData")
head(count_table)
is.matrix(count_table)
is.numeric(count_table)

pheno <- data.frame(tp = c("T1", "T1", "T1", "T24", "T24", "T24"))
deseq_dataset <- DESeqDataSetFromMatrix(countData = count_table ,
                                        colData = pheno , design = ~tp)
colData(deseq_dataset)

#Task 7.2
deseq_dataset <- estimateSizeFactors(deseq_dataset)
sizeFactors(deseq_dataset)
counts(deseq_dataset , normalized = TRUE)

#Task 7.3
deseq_dataset <- estimateDispersions(deseq_dataset)
plotDispEsts(deseq_dataset)

#Task 7.4
deseq_dataset <- nbinomWaldTest(deseq_dataset)
results_table <- results(deseq_dataset , contrast = c("tp", "T24","T1"))
head(results_table)
summary(results_table)

#Task 8
plotMA(deseq_dataset , alpha = 0.01)
plotDispEsts(deseq_dataset)

#log transformation
rlog_data <- rlogTransformation(deseq_dataset , blind = TRUE)

#Task 9

library('gplots')

dist_rl = dist(t(assay(rlog_data)))
dist_mat = as.matrix(dist_rl)
heatmap.2(dist_mat , trace = "none")

plotPCA(rlog_data , intgroup = "tp")

#Task 10
load(file = 'practical_2/annotation/annotation_table.RData')

dim(results_table)

results_table$PROKKA.ID <- row.names(results_table)
head(results_table)
res_table_2 = results_table[-c(0)]
head(res_table_2)
merge(as.data.frame(results_table), annotation_table, by="PROKKA.ID")

#Task 11
results_table$padj < 0.05 & abs(results_table$log2FoldChange) > 1
