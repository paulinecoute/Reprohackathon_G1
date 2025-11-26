# Calcul de corrélation entre notre table de comptage et celle des auteurs

# Notre table : data_counts/counts.txt 
# Table des auteurs : data_others/GSE139659_IPvsctrl.complete.txt
# Attention les tables n'ont pas la même structure (nombre et noms des colonnes etc)

# Chargement et Préparation de notre table 
us=read.table("data_counts/counts.txt", header = TRUE, sep = "\t", row.names = 1)
us=us[ , 6:ncol(us) ] # supression des colonnes Chr Start End Strand Length
rownames(us)=sub("gene-", "", rownames(us)) # enlever gene- au début des ID
colnames(us)=c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6") # renommage des colonnes pour coller aux auteurs 
# ATTENTION l'ordre est inversé par rapport aux auteurs, on a d'abord les IP puis les ctrl
nb_genes_us=nrow(us)

# Chargement et Préparation de la table des auteurs 
authors=read.table("data_others/GSE139659_IPvsctrl.complete.txt",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
authors=authors[, c("Name", "ctrl4", "ctrl5", "ctrl6", "IP1", "IP2", "IP3")] # garder que les colonnes nom + comptages
rownames(authors)=authors$Name
authors$Name=NULL
nb_genes_auth = nrow(authors)

# Mise en correspondance des tables 
common_genes=intersect(rownames(us), rownames(authors))
nb_common_genes=length(common_genes)
nb_common_genes
# Matrices 
us_common=us[common_genes, ]
auth_common=authors[common_genes, ]

# Calcul corrélation sur les gènes en commun 
cor_common=sapply(colnames(us_common), function(sample) {cor(us_common[, sample], auth_common[, sample], method = "pearson")}) 
# gènes en commun 
# une valeur par colonne 

# Résultats 
cat("Nombre de gènes dans notre table :", nb_genes_us, "\n") #2967
cat("Nombre de gènes dans la table auteurs :", nb_genes_auth, "\n") #2380
cat("Nombre de gènes en commun :", nb_common_genes, "\n") #2372
cat("Corrélation sur les gènes en commun :", cor_common, "\n") #0.9997867 0.9997522 0.9998411 0.9995568 0.9994292 0.9971949 

# Il est possible qu'on mette plutôt les résultats dans un rapport txt pour que ce soit bien intégré dans le pipeline
