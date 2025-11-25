# Fonctionne sur une machine locale 
# Calcul de corrélation entre notre table de comptage et celle des auteurs

# Notre table : data_counts/counts.txt 
# Table des auteurs : data_others/GSE139659_IPvsctrl.complete.txt
# Attention les tables n'ont pas la même structure (nombre et noms des colonnes etc)

################################################################################
# CHARGEMENT DES BIBLIOTHÈQUES
################################################################################

library(dplyr)



################################################################################
# PRÉTRAITEMENT DES DONNÉES
################################################################################

# Chargement et Préparation de notre table de comptage
us = read.table("../data_others/counts.txt", header = TRUE, sep = "\t", row.names = 1)
us = us[ , 6:ncol(us) ] # supression des colonnes Chr Start End Strand Length
rownames(us) = sub("gene-", "", rownames(us)) # enlever gene- au début des ID
colnames(us) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6") # renommage des colonnes pour coller aux auteurs 
# ATTENTION l'ordre est inversé par rapport aux auteurs, on a d'abord les IP puis les ctrl
us = us %>% select(4, 5, 6, 1, 2, 3)  # auteurs et us sont organisées de la même façon
nb_genes_us=nrow(us)


# table featurecount -s 0
s0 = read.table("../data_others/counts_s0.txt", header = TRUE, sep = "\t", row.names = 1)
s0 = s0[ , 6:ncol(s0) ] # supression des colonnes Chr Start End Strand Length
rownames(s0) = sub("gene-", "", rownames(s0)) # enlever gene- au début des ID
colnames(s0) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6") # renommage des colonnes pour coller aux auteurs 
# ATTENTION l'ordre est inversé par rapport aux auteurs, on a d'abord les IP puis les ctrl
s0 = s0 %>% select(4, 5, 6, 1, 2, 3)  # auteurs et us sont organisées de la même façon
nb_genes_s0=nrow(s0)


# table featurecount -s 1
s1 = read.table("../data_others/counts_s1.txt", header = TRUE, sep = "\t", row.names = 1)
s1 = s1[ , 6:ncol(s1) ] # supression des colonnes Chr Start End Strand Length
rownames(s1) = sub("gene-", "", rownames(s1)) # enlever gene- au début des ID
colnames(s1) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6") # renommage des colonnes pour coller aux auteurs 
# ATTENTION l'ordre est inversé par rapport aux auteurs, on a d'abord les IP puis les ctrl
s1 = s1 %>% select(4, 5, 6, 1, 2, 3)  # auteurs et us sont organisées de la même façon
nb_genes_s1=nrow(s1)


# table featurecount -s 2
s2 = read.table("../data_others/counts_s2.txt", header = TRUE, sep = "\t", row.names = 1)
s2 = s2[ , 6:ncol(s2) ] # supression des colonnes Chr Start End Strand Length
rownames(s2) = sub("gene-", "", rownames(s2)) # enlever gene- au début des ID
colnames(s2) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6") # renommage des colonnes pour coller aux auteurs 
# ATTENTION l'ordre est inversé par rapport aux auteurs, on a d'abord les IP puis les ctrl
s2 = s2 %>% select(4, 5, 6, 1, 2, 3)  # auteurs et us sont organisées de la même façon
nb_genes_s2=nrow(s2)


# Chargement et Préparation de la table des auteurs 
authors=read.csv("../data_others/GSE139659_IPvsctrl.complete.xls",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
authors=authors[, c("Name", "ctrl4", "ctrl5", "ctrl6", "IP1", "IP2", "IP3")] # garder que les colonnes nom + comptages
rownames(authors)=authors$Name
authors$Name=NULL
nb_genes_auth = nrow(authors)


# Chargement de la table des auteurs complete
complete =read.csv("../data_others/GSE139659_IPvsctrl.complete.xls",header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Mise en correspondance des tables 
common_genes=intersect(rownames(us), rownames(authors))
nb_common_genes=length(common_genes)
nb_common_genes


# Matrices 
# Pas nécessaire dans notre cas mais pour d'autres occasions oui
us_common=us[common_genes, ]
auth_common=authors[common_genes, ]

# featurecount
s0_common= s0[common_genes, ]
s1_common= s1[common_genes, ]
s2_common= s2[common_genes, ]

nrow(s0_common); nrow(s1_common); nrow(s2_common)
# 2967 gènes dans les tables s*_common donc mêmes gènes que dans la table authors



################################################################################
# ANALYSES DES DONNÉES - évaluation de la relation entre données authors & us 
################################################################################

# Calcul corrélation sur les gènes en commun 
cor_common=sapply(colnames(us_common), function(sample) {cor(us_common[, sample], auth_common[, sample], method = "pearson")})
# gènes en commun 
# gènes en commun 
# une valeur par colonne 

# featurecount
cor_s0=sapply(colnames(s0_common), function(sample) {cor(s0_common[, sample], auth_common[, sample], method = "pearson")})
cor_s1=sapply(colnames(s1_common), function(sample) {cor(s1_common[, sample], auth_common[, sample], method = "pearson")})
cor_s2=sapply(colnames(s2_common), function(sample) {cor(s2_common[, sample], auth_common[, sample], method = "pearson")})



################################################################################
# RÉSULTATS
################################################################################

cat("Nombre de gènes dans notre table :", nb_genes_us, "\n") #2967
cat("Nombre de gènes dans la table des auteurs :", nb_genes_auth, "\n") #2967
cat("Nombre de gènes en commun :", nb_common_genes, "\n") #2967
cat("Corrélation sur les gènes en commun : \n
     ctrl4     ctrl5     ctrl6       IP1       IP2       IP3 \n",
    cor_common, "\n") # 0.9996447 0.9997041 0.9974026 0.9997189 0.9997153 0.9996640

# feature count

cor_s0; cor_s1; cor_s2
