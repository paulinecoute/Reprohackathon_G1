# Reproduction des figures de l'article : 
# MA plot sur tous les gènes
# MA plot sur les gènes impliqués dans la traduction 

# Inputs 
# Table de comptage data_counts/counts.txt 
# Liste des gènes impliqués dans la traduction : data_others/genes_translation.txt

dir.create("plots") # dossier dans lequel seront mises les figures

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

dev.copy(png, filename = "plots/MAplot_all.png", width = 600, height = 500) # dimensions qui se rapporchent de celles des auteurs 
dev.off()

# MA plot des gènes impliqués dans la traduction 

# Chargement et Préparation du fichier avec les noms des gènes de la traduction 
translation=scan( "data_others/genes_translation.txt",what = "character", comment.char = "#", quiet = TRUE)
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

translation_lines=readLines("genes_translation.txt")
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

dev.copy(png, filename = "plots/MAplot_translation.png", width = 400, height = 400)
dev.off()



