##########################################################################

# Correlation analysis between our gene count table and the authors' dataset

# Inputs
# Our gene count table: data_counts/counts.txt
# Authors' gene count table: others/GSE139659_IPvsctrl.complete.xls 

# Outputs 
# Correlation report saved in the "reports" directory

##########################################################################

dir.create("reports")

# The two count tables use different formats, so we first harmonize them

us = read.table("data_counts/counts.txt", header = TRUE, sep = "\t", row.names = 1)
us = us[, 6:ncol(us)]           
rownames(us) = sub("gene-", "", rownames(us)) 
colnames(us) = c("IP1", "IP2", "IP3", "ctrl4", "ctrl5", "ctrl6")
sample_order = c("ctrl4", "ctrl5", "ctrl6", "IP1", "IP2", "IP3")
us = us[, sample_order]
nb_genes_us = nrow(us)

authors = read.csv("others/GSE139659_IPvsctrl.complete.xls",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
authors = authors[, c("Name", sample_order)] 
rownames(authors) = authors$Name
authors$Name = NULL
authors = authors[, sample_order]            
nb_genes_auth = nrow(authors)

# Correlation

common_genes = intersect(rownames(us), rownames(authors))
nb_common_genes = length(common_genes)
us_common = us[common_genes, ]
auth_common = authors[common_genes, ]
cor_common = sapply(sample_order, function(sample) {cor(us_common[, sample], auth_common[, sample], method = "pearson")})

sink("reports/correlation.txt")
cat("Nombre de gènes dans notre table :", nb_genes_us, "\n")
cat("Nombre de gènes dans la table auteurs :", nb_genes_auth, "\n")
cat("Nombre de gènes en commun :", nb_common_genes, "\n")
cat("Corrélation sur les gènes en commun :\n")
print(cor_common)
sink()

