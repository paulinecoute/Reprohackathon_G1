
# script to perform differential gene expression analysis using DESeq2 package
# Les gènes ont été comptés à l'aide de featureCounts version 1.4.6-p3 du package Subreads (paramètres : -t gene -g ID -s 1). Les données de comptage ont été analysées à l'aide de R version 3.4.1 (doi : 10.1093/bioinformatics/btt656) avec le package Bioconductor DESeq2 version 1.16 et le package EnrichmentBrowser Rversion 2.14.3 [https://www.R-project.org/].

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: preparing count data ----------------
data(airway)

sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)


# read in sample info
colData <- read.csv('sample_info.csv')


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated

# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# MA plot
plotMA(res)
