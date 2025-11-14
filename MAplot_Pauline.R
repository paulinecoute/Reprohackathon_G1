library(DESeq2)
library(tidyverse)

#############################################
### 1) Charger ton fichier counts annotés
#############################################

counts <- read_delim("counts_annotated.txt", delim="\t")

# Extraire la matrice de comptage (seulement les colonnes des 6 conditions)
count_matrix <- counts %>%
  select(ctrl4, ctrl5, ctrl6, IP1, IP2, IP3)

# Geneid comme rownames
rownames(count_matrix) <- counts$Geneid

#############################################
### 2) Construire colData (design)
#############################################

colData <- data.frame(
  condition = factor(c("ctrl","ctrl","ctrl","IP","IP","IP"))
)
rownames(colData) <- colnames(count_matrix)

#############################################
### 3) Construire l'objet DESeq2
#############################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)

#############################################
### 4) Pré-filtrage
#############################################
dds <- dds[rowSums(counts(dds)) >= 10, ]

#############################################
### 5) Définir la condition de référence
#############################################
dds$condition <- relevel(dds$condition, ref="ctrl")

#############################################
### 6) DESeq
#############################################
dds <- DESeq(dds)
res <- results(dds)

#############################################
### 7) MA plot
#############################################
plotMA(res, ylim=c(-4,4))
