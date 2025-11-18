library(DESeq2)
library(tidyverse)

counts <- read_delim("counts_annotated.txt", delim="\t")

count_matrix <- counts %>%
  select(ctrl4, ctrl5, ctrl6, IP1, IP2, IP3) # que les colonnes des comptages 

rownames(count_matrix) <- counts$Geneid

colData <- data.frame(condition = factor(c("ctrl","ctrl","ctrl","IP","IP","IP")))
rownames(colData) <- colnames(count_matrix)

# objet deseq2
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ] #prefiltrer 

dds$condition <- relevel(dds$condition, ref="ctrl") # controle est la ref f


dds <- DESeq(dds)
res <- results(dds)

#ma plot 
plotMA(res, ylim=c(-4,4))
