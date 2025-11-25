rm(list=objects());graphics.off()

# =========================
# Inputs 
# =========================
# - Table de comptage: data_counts/counts.txt 
# - Liste des gènes impliqués dans la traduction : data_others/genes_translation.txt

dir.create("plots", showWarnings = FALSE)

# =========================
# Load Count Table
# =========================
count_data = read.table("data_others/counts.txt", header = TRUE, row.names = 1)

# Keep only count columns
count_data = count_data[ , 6:ncol(count_data) ]

# Metadata
metadata = data.frame(condition = factor(c(rep("Control", 3), rep("Treated", 3))))
rownames(metadata) = colnames(count_data)
metadata$condition = relevel(metadata$condition, ref = "Treated")

# =========================
# DESeq2
# =========================
library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
dds = DESeq(dds)
res = results(dds)
res$GeneID = sub("^gene-", "", rownames(res))

# =========================
# VOLCANO PLOT #1 – ALL GENES
# =========================

res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# SIGNIFICANT GENES
sig <- res_df[res_df$padj < 0.05, ]

up  <- sum(sig$log2FoldChange > 0)
down <- sum(sig$log2FoldChange < 0)

cat("Upregulated genes:", up, "\n")
cat("Downregulated genes:", down, "\n")


png("plots/Volcano_all.png", width = 650, height = 500)

with(res, plot(
  x = log2FoldChange,
  y = -log10(padj),
  pch = 20,
  cex = 0.6,
  col = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "grey"),
  xlab = "Log2 Fold Change (Treated vs Control)",
  ylab = "-log10(FDR)",
  xlim = c(-6, 6),
  main = paste0(
    "Volcano plot of complete RNA-seq dataset\n",
    "Up: ", up, "   Down: ", down
  )
))

abline(v = c(-1, 1), col = "darkgrey", lty = 2)
abline(h = -log10(0.05), col = "darkgrey", lty = 2)

dev.off()


# =========================
# Load Translation Gene List
# =========================
translation = scan("data_others/genes_translation.txt", what = "character",
                   comment.char = "#", quiet = TRUE)

translation_ids = translation[grepl("^SAOUHSC_", translation)]

# Extract DESeq2 results for translation genes
res_translation = res[res$GeneID %in% translation_ids, ]
cat("Translation genes found in dataset:", nrow(res_translation), "\n")


# =========================
# Parse annotation file for labeling
# =========================
translation_lines = readLines("data_others/genes_translation.txt")
translation_data = translation_lines[grepl("^SAOUHSC_", translation_lines)]

parsed_translation <- do.call(rbind, lapply(translation_data, function(line) {
  parts <- strsplit(line, " ")[[1]]
  locus_tag <- parts[1]
  symbol <- gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

target_symbols <- c("frr", "infA", "tsf", "infC", "infB", "pth")

genes_to_label = parsed_translation[parsed_translation$symbol %in% target_symbols, ]
missing = setdiff(target_symbols, genes_to_label$symbol)

if (length(missing) > 0) {
  missing_df = data.frame(
    locus_tag = NA,
    symbol = missing,
    stringsAsFactors = FALSE
  )
  genes_to_label = rbind(genes_to_label, missing_df)
}

annot = merge(as.data.frame(res_translation), genes_to_label,
              by.x = "GeneID", by.y = "locus_tag", all.y = TRUE)



# =========================
# Identify AA–tRNA synthetases
# =========================
aa_trna_lines = translation_lines[grepl("tRNA synthetase", translation_lines)]

aa_trna_data = do.call(rbind, lapply(aa_trna_lines, function(line) {
  parts = strsplit(line, " ")[[1]]
  locus_tag = parts[1]
  symbol = gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

res_translation$AA_tRNA = res_translation$GeneID %in% aa_trna_data$locus_tag

# =========================
# Top 10 most significant translation genes
# =========================

res_translation_df <- as.data.frame(res_translation)
res_translation_df <- res_translation_df[!is.na(res_translation_df$padj), ]

# Select 10 smallest FDR values
top10 <- res_translation_df[order(res_translation_df$padj), ][1:10, ]

# Prepare table for labeling
top10$label = top10$GeneID
# If symbol exists in parsed data, use it
top10 = merge(top10, parsed_translation, by.x="GeneID", by.y="locus_tag", all.x=TRUE)
top10$label = ifelse(is.na(top10$symbol), top10$label, top10$symbol)


# =========================
# VOLCANO PLOT #2 – TRANSLATION GENES
# =========================
png("plots/Volcano_translation.png", width = 500, height = 500)

with(res_translation, plot(
  x = log2FoldChange,
  y = -log10(padj),
  pch = 20,
  cex = 0.7,
  col = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "red", "grey"),
  xlab = "Log2 Fold Change (Treated vs Control)",
  ylab = "-log10(FDR)",
  main = "Volcano Plot – Translation Genes",
  xlim = c(-6, 6),
  ylim = c(0, max(-log10(res_translation$padj), na.rm = TRUE) + 1)
))

abline(v = c(-1, 1), col = "darkgrey", lty = 2)
abline(h = -log10(0.05), col = "darkgrey", lty = 2)

# Add labels for top 10 most significant genes
for (i in 1:nrow(top10)) {
  x = top10$log2FoldChange[i]
  y = -log10(top10$padj[i])
  text(x, y + 0.5, top10$label[i], cex = 0.8, font = 2, col="green")
  arrows(x0 = x, y0 = y + 0.4, x1 = x, y1 = y, length = 0.08)
}

# Add labels for target genes
for (i in 1:nrow(annot)) {
  if (!is.na(annot$GeneID[i])) {
    x = annot$log2FoldChange[i]
    y = -log10(annot$padj[i])
    text(x, y + 0.4, annot$symbol[i], cex = 0.85, font = 2)
    arrows(x0 = x, y0 = y + 0.3, x1 = x, y1 = y, length = 0.08)
  } else {
    message("Gene not found: ", annot$symbol[i])
  }
}

# Highlight AA–tRNA synthetases
points(
  x = res_translation$log2FoldChange[res_translation$AA_tRNA],
  y = -log10(res_translation$padj[res_translation$AA_tRNA]),
  pch = 21,
  bg = ifelse(res_translation$padj[res_translation$AA_tRNA] < 0.05, "red", "grey"),
  col = "black",
  lwd = 1.5,
  cex = 1
)

legend("topright", legend = "AA-tRNA synthetase", pch = 21,
       pt.bg = "grey", col = "black", bty = "n", cex = 0.8)

dev.off()



