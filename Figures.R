##########################################################################

# Differential expression analysis from the gene count table using DESeq2
# Computation of the number of up- and down-regulated genes
# Generation of figures: MA plots and volcano plots

# Inputs
# Gene count table: data_counts/counts.txt
# List of genes involved in translation: others/genes_translation.txt

# Outputs 
# All figures are saved in the "figures" directory
# Differential analysis report saved in "reports"

##########################################################################

dir.create("figures") 

# Load and format the gene count table

count_data=read.table("data_counts/counts.txt", header = TRUE, row.names = 1)
count_data=count_data[ , 6:ncol(count_data) ] 
metadata=data.frame(condition = factor(c(rep("Treated", 3), rep("Control", 3)))) 
rownames(metadata)=colnames(count_data)
metadata$condition = relevel(metadata$condition, ref = "Control") 

# Differential Analysis with DESeq2 

library(DESeq2)
dds=DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
dds=DESeq(dds)
res=results(dds)
genes_up = rownames(res[ !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 0 , ])
genes_down = rownames(res[ !is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < 0 , ])
genes_ns = rownames(res[ is.na(res$padj) | res$padj >= 0.05 , ])

sink("reports/diff_analysis.txt")
cat("Number of UP-regulated genes: ", length(genes_up), "\n")
cat("Number of DOWN-regulated genes: ", length(genes_down), "\n")
cat("Number of non-significant genes: ", length(genes_ns), "\n\n")
cat("List of UP-regulated genes:\n")
print(genes_up)
cat("\nList of DOWN-regulated genes:\n")
print(genes_down)
cat("\nList of non-significant genes:\n")
print(genes_ns)
sink()

# MA plots: all genes and genes involved in translation

png("figures/MAplot_all.png", width = 600, height = 500)

plot(
  x = res$baseMean,
  y = res$log2FoldChange,
  log = "x",
  pch = 20,
  cex = 0.5,
  col = ifelse(res$padj < 0.05, "red", "black"),
  xlab = "Mean of normalized counts",
  ylab = "Log2 fold change",
  main = "MA-plot of complete RNA-seq dataset",
  ylim = c(-4, 4),
  xaxt = "n"
)

axis(1, at = c(1, 100, 10000), labels = expression(10^0, 10^2, 10^4))
abline(h = 0, col = "black", lty = 2)

dev.off()

# Load and format the gene involved in translation table
translation=scan( "others/genes_translation.txt",what = "character", comment.char = "#", quiet = TRUE)
translation_ids = translation[grepl("^SAOUHSC_", translation)] 
res$GeneID = sub("^gene-", "", rownames(res))  
res_translation = res[res$GeneID %in% translation_ids, ] 
genestranslation = res_translation$GeneID 
log2_baseMean=log2(res_translation$baseMean + 1) 
translation_up = res_translation$GeneID[ res_translation$padj < 0.05 & res_translation$log2FoldChange > 0 ]
translation_down = res_translation$GeneID[ res_translation$padj < 0.05 & res_translation$log2FoldChange < 0 ]
translation_ns = res_translation$GeneID[ is.na(res_translation$padj) | res_translation$padj >= 0.05 ]

sink("reports/translation.txt")
cat("Total translation-related genes detected: ", nrow(res_translation), "\n")
cat("Number of UP-regulated translation genes: ", length(translation_up), "\n")
cat("Number of DOWN-regulated translation genes: ", length(translation_down), "\n")
cat("Number of non-significant translation genes: ", length(translation_ns), "\n\n")
cat("List of ALL translation-related genes detected:\n")
print(genestranslation)
cat("\nList of UP-regulated translation genes:\n")
print(translation_up)
cat("\nList of DOWN-regulated translation genes:\n")
print(translation_down)
cat("\nList of non-significant translation genes:\n")
print(translation_ns)
sink()

png("figures/MAplot_translation.png", width = 350, height = 400)

plot(
  x = log2_baseMean,
  y = res_translation$log2FoldChange,
  pch = 20,
  cex = 0.7,
  col = ifelse(res_translation$padj < 0.05, "red", "grey"),
  xlab = "Log2 base Mean",
  ylab = "Log2 Fold Change",
  main = "MA-plot of genes related to translation",
  xlim = c(0, 20),
  ylim = c(-6, 5),
  xaxt = "n",
  yaxt = "n"
)

axis(1, at = seq(0, 20, by = 2))
axis(2, at = seq(-6, 5, by = 1))
abline(h = 0, col = "black", lty = 2)
legend("bottomleft",legend = c("Significant", "Non-significant"),col = c("red", "grey"),pch = 20,bty = "n",cex = 0.8)

translation_lines=readLines("others/genes_translation.txt")
translation_data=translation_lines[grepl("^SAOUHSC_", translation_lines)]

parsed_translation = do.call(rbind, lapply(translation_data, function(line) {
  parts = strsplit(line, " ")[[1]]
  locus_tag = parts[1]
  symbol = gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

target_symbols = c("frr", "infA", "tsf", "infC", "infB", "pth") # genes of interest (from the reference figure)
genes_to_label = parsed_translation[parsed_translation$symbol %in% target_symbols, ]
genes_to_label = rbind(genes_to_label,data.frame(locus_tag = "SAOUHSC_00475",symbol = "pth",stringsAsFactors = FALSE))

missing=setdiff(target_symbols, genes_to_label$symbol)
if (length(missing) > 0) {
  missing_df=data.frame(locus_tag = NA,symbol = missing,stringsAsFactors = FALSE)
  genes_to_label=rbind(genes_to_label, missing_df)
}

annot=merge(as.data.frame(res_translation),genes_to_label,by.x = "GeneID",by.y = "locus_tag",all.y = TRUE)

offsets = list(
  pth  = c(dx = -1.6, dy = -1.6),   
  frr  = c(dx = -1.6, dy =  1.6),   
  infA = c(dx = -1.6, dy =  1.6),   
  tsf  = c(dx =  1.6, dy =  1.6),   
  infC = c(dx =  2.0, dy = -0.8),
  infB = c(dx =  1.6, dy = -1.6)    
)

for (i in 1:nrow(annot)) {
  
  if (!is.na(annot$GeneID[i])) {
    
    x = log2(annot$baseMean[i] + 1)
    y = annot$log2FoldChange[i]
    gene = annot$symbol[i]
    
    if (!(gene %in% names(offsets))) next
    
    dx = offsets[[gene]]["dx"]
    dy = offsets[[gene]]["dy"]

    xf = x + dx
    yf = y + dy
    xt = x + dx * 1.25
    yt = y + dy * 1.25
    
    segments(x0 = x,  y0 = y,x1 = xf, y1 = yf,length = 0.10,angle = 20,lwd = 1.5)
    
    text(xt, yt,labels = gene,cex = 1, font = 2)
  }
}

aa_trna_lines = translation_lines[grepl("tRNA synthetase", translation_lines)]

aa_trna_data = do.call(rbind, lapply(aa_trna_lines, function(line) {
  parts = strsplit(line, " ")[[1]]
  locus_tag = parts[1]
  symbol = gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

res_translation$AA_tRNA = res_translation$GeneID %in% aa_trna_data$locus_tag

points(
  x = log2(res_translation$baseMean[res_translation$AA_tRNA] + 1),
  y = res_translation$log2FoldChange[res_translation$AA_tRNA],
  pch = 21,
  bg = ifelse(res_translation$padj[res_translation$AA_tRNA] < 0.05, "red", "grey"),
  col = "black",
  lwd = 1.5,
  cex = 0.9
)

usr=par("usr") 
legend(x = usr[2] - 10,  y = usr[3] + 1,legend = c("AA-tRNA synthetase"),pch = 21,pt.bg = "grey",col = "black",lwd = 1.5,pt.cex = 1,bty = "n",cex = 0.8 )

dev.off()

# Volcano plots: all genes and genes involved in translation 

res_df = as.data.frame(res)
res_df = res_df[!is.na(res_df$padj), ]

sig = res_df[res_df$padj < 0.05, ]
up = sum(sig$log2FoldChange > 0)
down = sum(sig$log2FoldChange < 0)

png("figures/Volcano_all.png", width = 650, height = 500)

with(res_df, plot(
  x = log2FoldChange,
  y = -log10(padj),
  pch = 20,
  cex = 0.6,
  col = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "grey"),
  xlab = "Log2 Fold Change (Treated vs Control)",
  ylab = "-log10(FDR)",
  xlim = c(-6, 6),
  main = paste0("Volcano plot – All genes\nUp: ", up, "   Down: ", down)
))

abline(v = c(-1, 1), col = "darkgrey", lty = 2)
abline(h = -log10(0.05), col = "darkgrey", lty = 2)

dev.off()

res_translation_df = as.data.frame(res_translation)
res_translation_df = res_translation_df[!is.na(res_translation_df$padj), ]

# Top 10
top10 = res_translation_df[order(res_translation_df$padj), ][1:10, ]
top10 = merge(top10, parsed_translation, by.x = "GeneID", by.y = "locus_tag", all.x = TRUE)
top10$label = ifelse(is.na(top10$symbol), top10$GeneID, top10$symbol)

png("figures/Volcano_translation.png", width = 500, height = 500)

with(res_translation_df, plot(
  x = log2FoldChange,
  y = -log10(padj),
  pch = 20,
  cex = 0.7,
  col = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "grey"),
  xlab = "Log2 Fold Change (Treated vs Control)",
  ylab = "-log10(FDR)",
  main = "Volcano plot – Translation genes",
  xlim = c(-6, 6),
  ylim = c(0, max(-log10(padj), na.rm = TRUE) + 1)
))

abline(v = c(-1, 1), col = "darkgrey", lty = 2)
abline(h = -log10(0.05), col = "darkgrey", lty = 2)

for (i in 1:nrow(top10)) {
  x = top10$log2FoldChange[i]
  y = -log10(top10$padj[i])
  text(x, y + 0.5, top10$label[i], cex = 0.8, font = 2, col="green")
  arrows(x0 = x, y0 = y + 0.4, x1 = x, y1 = y, length = 0.08)
}

for (i in 1:nrow(genes_to_label)) {
  gene_id <- genes_to_label$locus_tag[i]
  row <- res_translation_df[res_translation_df$GeneID == gene_id, ]
  if (nrow(row) == 1) {
    x = row$log2FoldChange
    y = -log10(row$padj)
    text(x, y + 0.4, genes_to_label$symbol[i], cex = 0.85, font = 2)
    arrows(x0 = x, y0 = y + 0.3, x1 = x, y1 = y, length = 0.08)
  }
}

points(
  x = res_translation_df$log2FoldChange[res_translation_df$AA_tRNA],
  y = -log10(res_translation_df$padj[res_translation_df$AA_tRNA]),
  pch = 21,
  bg = ifelse(res_translation_df$padj[res_translation_df$AA_tRNA] < 0.05, "red", "grey"),
  col = "black",
  lwd = 1.5,
  cex = 1
)

legend("topright", legend = "AA-tRNA synthetase", pch = 21,
       pt.bg = "grey", col = "black", bty = "n", cex = 0.8)

dev.off()

# Heatmap and boxplots for regulated translation genes

library(pheatmap)
library(reshape2)
library(ggplot2)


norm_counts = counts(dds, normalized = TRUE)
norm_log = log2(norm_counts + 1)
sig_translation = res_translation[!is.na(res_translation$padj) & res_translation$padj < 0.05, ]
genes_sig_ids = sig_translation$GeneID
genes_sig_full = paste0("gene-", genes_sig_ids)
genes_sig_full = genes_sig_full[genes_sig_full %in% rownames(norm_log)]
mat = norm_log[genes_sig_full, , drop = FALSE]

png("figures/Heatmap_translation_sig.png", width = 700, height = 600)
pheatmap(
  mat,
  annotation_col = metadata,
  scale = "row",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  main = "Heatmap – Significant translation genes"
)
dev.off()


genes_up_ids   = sig_translation$GeneID[sig_translation$log2FoldChange > 0]
genes_down_ids = sig_translation$GeneID[sig_translation$log2FoldChange < 0]
genes_up_full   = paste0("gene-", genes_up_ids)
genes_down_full = paste0("gene-", genes_down_ids)
genes_up_full   = intersect(genes_up_full, rownames(norm_log))
genes_down_full = intersect(genes_down_full, rownames(norm_log))


get_df = function(gene_list, direction_label) {
  if (length(gene_list) == 0) return(NULL)
  
  df = as.data.frame(t(norm_log[gene_list, , drop = FALSE]))
  df$condition = metadata$condition
  
  df_m = melt(df, id.vars = "condition", variable.name = "gene", value.name = "expression")
  df_m$direction = direction_label
  return(df_m)
}

df_up = get_df(genes_up_full, "Up-regulated")
df_down = get_df(genes_down_full, "Down-regulated")
df_all = rbind(df_up, df_down)

df_all$group = paste(df_all$direction, df_all$condition, sep = "_")

# Order groups visually
df_all$group = factor(df_all$group,levels = c("Up-regulated_Control","Up-regulated_Treated","Down-regulated_Control","Down-regulated_Treated"))

png("figures/Boxplot_translation_sig.png", width = 650, height = 500)

ggplot(df_all, aes(x = group, y = expression, fill = direction)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c("Up-regulated" = "#1f77b4", 
                               "Down-regulated" = "#ff7f0e")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  xlab("") +
  ylab("Log2 normalized expression") +
  ggtitle("Expression of translation genes – Up/Down × Control/Treated")

dev.off()




