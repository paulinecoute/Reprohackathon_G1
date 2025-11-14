# Script R qui ajoute dans counts.txt les annotations de reference.gff pour créer un nouveau tableau et construire après les MA plot 
# Pour l'instant il s'appelle seul ce n'est pas automatisé donc on le fera plus tard une fois que tout sera validé
# Dernère mise-à-jour : 14.11.25 


# Packages nécessaires 
library(dplyr)
library(readr)
library(stringr)

# Table reference.gff
gff <- read_delim("reference.gff",delim = "\t",comment = "#",col_names = FALSE)
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
cds <- gff %>% filter(type == "CDS") # on garde que les CDS 

# On extrait locus_tag, protein_id, product (en vrai c'est surtout protein_id et product qui nous intéresse mais je me suis dit qu'avoir le locus tag ça peut être pas mal aussi)

extract_attr <- function(attr, key) {
  pattern <- paste0(key, "=([^;]+)")
  str_match(attr, pattern)[,2]
}

annot <- cds %>%
  mutate(
    locus_tag  = extract_attr(attributes, "locus_tag"),
    protein_id = extract_attr(attributes, "protein_id"),
    product    = extract_attr(attributes, "product")
  ) %>%
  select(locus_tag, protein_id, product) %>%
  distinct() %>%
  filter(!is.na(locus_tag))

# Table counts.txt 
counts <- read_delim("counts.txt",delim = "\t",comment = "#")
counts$Geneid <- str_replace(counts$Geneid, "^gene-", "") # on enlève juste gene- qu'il y a au début 

# Fusion des 2 tables 
counts_annot <- counts %>%
  left_join(annot, by = c("Geneid" = "locus_tag"))

# On ne conserve que les ID, product (pour plus tard filtrer les gènes impliqués dans la traduction) et les comptages 
counts_clean <- counts_annot %>%
  select(Geneid, protein_id, product, matches("SRR"))

# Renommage des colonnes pour se rapprocher de la table des auteurs 
new_colnames <- c("ctrl4", "ctrl5", "ctrl6", "IP1", "IP2", "IP3")
colnames(counts_clean)[4:9] <- new_colnames

# Nouveau tableau annoté 
write_delim(counts_clean,"counts_annotated.txt",delim = "\t")
