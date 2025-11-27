# Notre table : data_counts/counts.txt 
# Table des auteurs : others/GSE139659_IPvsctrl.complete.xls (TSV renommé)

dir.create("correlation")

# Chargement et Préparation de notre table 
us = read.table("data_counts/counts.txt",
                header = TRUE, sep = "\t", row.names = 1)
us = us[, 6:ncol(us)]           # suppr Chr Start End Strand Length
rownames(us) = sub("gene-", "", rownames(us)) 
colnames(us) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6")

# ordre qu'on va utiliser PARTOUT (comme les auteurs)
sample_order = c("ctrl4", "ctrl5", "ctrl6", "IP1", "IP2", "IP3")
us = us[, sample_order]
nb_genes_us = nrow(us)

# Chargement et Préparation de la table des auteurs (TSV même si .xls)
authors = read.csv("others/GSE139659_IPvsctrl.complete.xls",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)


authors = authors[, c("Name", sample_order)]  # Name + mêmes colonnes
rownames(authors) = authors$Name
authors$Name = NULL
authors = authors[, sample_order]             # on force aussi l'ordre ici
nb_genes_auth = nrow(authors)

# Mise en correspondance des tables 
common_genes = intersect(rownames(us), rownames(authors))
nb_common_genes = length(common_genes)

# Matrices 
us_common   = us[common_genes, ]
auth_common = authors[common_genes, ]

# Calcul corrélation sur les gènes en commun 
cor_common = sapply(sample_order, function(sample) {
  cor(us_common[, sample], auth_common[, sample], method = "pearson")
})

# Résultats 
sink("correlation/correlation.txt")
cat("Nombre de gènes dans notre table :", nb_genes_us, "\n")
cat("Nombre de gènes dans la table auteurs :", nb_genes_auth, "\n")
cat("Nombre de gènes en commun :", nb_common_genes, "\n")
cat("Corrélation sur les gènes en commun :\n")
print(cor_common)
sink()

