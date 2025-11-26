# Reproduction des figures de l'article : 
# MA plot sur tous les gènes
# MA plot sur les gènes impliqués dans la traduction 

# un volcano plot est également effectué 

# Inputs 
# Table de comptage data_counts/counts.txt 
# Liste des gènes impliqués dans la traduction : others/genes_translation.txt

dir.create("figures") # dossier dans lequel seront mises les figures

# Chargement et Préparation de notre table 
count_data=read.table("data_counts/counts.txt", header = TRUE, row.names = 1)

# Garder seulement les colonnes de comptage (pas Chr, Start, etc.)
count_data=count_data[ , 6:ncol(count_data) ] # supression des colonnes Chr Start End Strand Length
metadata=data.frame(condition = factor(c(rep("Control", 3), rep("Treated", 3)))) # ajout des conditions
rownames(metadata)=colnames(count_data)
metadata$condition = relevel(metadata$condition, ref = "Treated") # fixation de la référence

# DESeq2 
library(DESeq2)
dds=DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ condition)
dds=DESeq(dds)
res=results(dds)

# MA plot de l'ensemble des gènes 

plot(
  x = res$baseMean,
  y = res$log2FoldChange,
  log = "x",
  pch = 20,
  cex = 0.5,
  col = ifelse(res$padj < 0.05, "red", "black"),
  xlab = "Mean of normalized counts", #log scale
  ylab = "Log2 fold change",
  main = "MA-plot of complete RNA-seq dataset",
  ylim = c(-4, 4),
  xaxt = "n"  # empêche le tracé automatique de l’axe x
)

# Mêmes axes / ligne horizontale que les auteurs 
axis(1, at = c(1, 100, 10000), labels = expression(10^0, 10^2, 10^4))
abline(h = 0, col = "black", lty = 2)

dev.copy(png, filename = "figures/MAplot_all.png", width = 600, height = 500) # dimensions qui se rapporchent de celles des auteurs 
dev.off()

# MA plot des gènes impliqués dans la traduction 

# Chargement et Préparation du fichier avec les noms des gènes de la traduction 
translation=scan( "others/genes_translation.txt",what = "character", comment.char = "#", quiet = TRUE)
translation_ids = translation[grepl("^SAOUHSC_", translation)] # valides, sans les com etc c'est juste une vérification supplémentaire 
res$GeneID = sub("^gene-", "", rownames(res))  # enlever "gene-" de l'ID
res_translation = res[res$GeneID %in% translation_ids, ] # on sélectionne les gènes présents dans notre table de comptage
genestranslation <- res_translation$GeneID # liste des gènes de la traduction
cat("Nombre de gènes impliqués dans la traduction présents dans notre table de comptage : ", nrow(res_translation), "\n") #172
log2_baseMean=log2(res_translation$baseMean + 1) # baseMean transformé en log2

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

# Mêmes axes / ligne horizontale / légende que les auteurs 
axis(1, at = seq(0, 20, by = 2))
axis(2, at = seq(-6, 5, by = 1))
abline(h = 0, col = "black", lty = 2)
legend("bottomleft",legend = c("Significant", "Non-significant"),col = c("red", "grey"),pch = 20,bty = "n",cex = 0.8)

# Ajout du nom des quelques gènes présents sur le graphique des auteurs 

translation_lines=readLines("others/genes_translation.txt")
translation_data=translation_lines[grepl("^SAOUHSC_", translation_lines)]

parsed_translation <- do.call(rbind, lapply(translation_data, function(line) {
  parts <- strsplit(line, " ")[[1]]
  locus_tag <- parts[1]
  symbol <- gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

target_symbols <- c("frr", "infA", "tsf", "infC", "infB", "pth") # gènes d'intérêts, présents sur la figure de l'article 
genes_to_label = parsed_translation[parsed_translation$symbol %in% target_symbols, ]

# NA si non trouvé 
missing=setdiff(target_symbols, genes_to_label$symbol)
if (length(missing) > 0) {
  missing_df=data.frame(
    locus_tag = NA,
    symbol = missing,
    stringsAsFactors = FALSE
  )
  genes_to_label=rbind(genes_to_label, missing_df)
}

# fusiion DESeq2 et annotation sur le graph
annot=merge(as.data.frame(res_translation),genes_to_label,by.x = "GeneID",by.y = "locus_tag",all.y = TRUE)

for (i in 1:nrow(annot)) {
  if (!is.na(annot$GeneID[i])) {
    x <- log2(annot$baseMean[i] + 1)
    y <- annot$log2FoldChange[i]
    
    text(x, y + 1.2, annot$symbol[i], cex = 0.9, font = 2)
    arrows(x0 = x, y0 = y + 1.1, x1 = x, y1 = y, length = 0.08)
  } else {
    message("Gène non trouvé : ", annot$symbol[i])
  }
}

# Ajout des AA-tRNA synthetases sur le graphe

aa_trna_lines = translation_lines[grepl("tRNA synthetase", translation_lines)]

aa_trna_data = do.call(rbind, lapply(aa_trna_lines, function(line) {
  parts = strsplit(line, " ")[[1]]
  locus_tag = parts[1]
  symbol = gsub(";", "", parts[2])
  data.frame(locus_tag = locus_tag, symbol = symbol, stringsAsFactors = FALSE)
}))

res_translation$AA_tRNA = res_translation$GeneID %in% aa_trna_data$locus_tag

# on ajoute une bordure noire sur ces points 
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

dev.copy(png, filename = "figures/MAplot_translation.png", width = 400, height = 400)
dev.off()

## Volcano plots 

# Dataset propre
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# Comptage up/down
sig <- res_df[res_df$padj < 0.05, ]
up  <- sum(sig$log2FoldChange > 0)
down <- sum(sig$log2FoldChange < 0)

##############
# Volcano all
##############
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

##############################
# Prépa volcano translation
##############################

res_translation_df <- as.data.frame(res_translation)
res_translation_df <- res_translation_df[!is.na(res_translation_df$padj), ]

# Top 10
top10 <- res_translation_df[order(res_translation_df$padj), ][1:10, ]
top10 <- merge(top10, parsed_translation, by.x = "GeneID", by.y = "locus_tag", all.x = TRUE)
top10$label <- ifelse(is.na(top10$symbol), top10$GeneID, top10$symbol)

###########################
# Volcano translation
###########################
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

# Labels top10
for (i in 1:nrow(top10)) {
  x = top10$log2FoldChange[i]
  y = -log10(top10$padj[i])
  text(x, y + 0.5, top10$label[i], cex = 0.8, font = 2, col="green")
  arrows(x0 = x, y0 = y + 0.4, x1 = x, y1 = y, length = 0.08)
}

# Labels gènes d’intérêt (target_symbols)
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

# Points aa–tRNA synthetases
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

