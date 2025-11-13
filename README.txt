Projet Reprohackathon M2 AMI2B
Pauline Couté, Ange-Louis Sammarcelli, Miguel Senovilla Herrero
Dernière mise-à-jour : 13.11.25 


DESCRIPTION DU PROJET 

# ajouter description brève du projet, but avec l’article etc 


ENVIRONNEMENT
 Le pipeline a été développé sur une machine virtuelle Ubuntu 24.04 du cloud IFB.
Pour garantir la compatibilité, il est fortement recommandé d’utiliser :
* Ubuntu 24.04
* Au minimum : 4 CPU
* Espace disque recommandé : 100 Go


INSTALLATIONS REQUISES 

Le workflow utilise Snakemake et Apptainer pour exécuter les outils dans des conteneurs générés à partir d’images Docker hébergées sur Docker Hub.

1. Copier le projet depuis votre machine locale vers la VM (Remplacez l’adresse IP par celle de votre VM) : scp -r ~/Reprohackathon_G1 ubuntu@11.22.33.44:~

2. Se connecter à la VM et mettre à jour les paquets : sudo apt update -y

3. Installer Snakemake sudo apt install -y snakemake

4. Installer Apptainer  sudo apt install -y build-essential libseccomp-dev pkg-config squashfs-tools cryptsetup uidmap
wget https://github.com/apptainer/apptainer/releases/download/v1.4.4/apptainer_1.4.4_amd64.deb
sudo dpkg -i apptainer_1.4.4_amd64.deb
sudo apt --fix-broken install -y

A titre informatif, les versions de Snakemake et Apptainer utilisées lors du développement sont respectivement 7.32.4 et apptainer version 1.4.4.  

LANCEMENT DU SCRIPT  
Se placer dans le dossier du projet et lancer la commande suivante :
snakemake --use-singularity -j 4

Le temps d’exécution du pipeline entier est d’environ 2h.

