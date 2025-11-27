Authors: Pauline Couté, Ange-Louis Sammarcelli, Miguel Senovilla Herrero
Last update: 13 November 2025


----- Project Description 

The goal of this project was to reproduce part of the analysis from Peyrusson et al., Nature Communications (2020), specifically, the RNA-seq workflow and several key figures presented by the authors. The entire analysis was implemented using Snakemake as the workflow management system and Docker for containerisation.

On the one hand, the pipeline performs the basic steps of an RNA-seq analysis, resulting in a gene count table. On the other hand, it carries out differential expression analysis using the DESeq2 R package and generates figures as well as statistics to compare our results with those reported by Peyrusson et al..


----- Environment 

The pipeline was developed on an Ubuntu 24.04 virtual machine provided by the IFB Cloud. To ensure full compatibility, we recommend running the workflow on Ubuntu 24.04 (or an equivalent Linux system), with at least 4 CPU cores and a minimum of 100 GB of available disk space.


----- Required Tools

The workflow relies on Snakemake and Apptainer to run all tools inside containers generated from Docker images hosted on Docker Hub.
To reproduce the environment, the project must first be transferred to the virtual machine. From your local machine, you can copy the repository as follows:
scp -r ~/Reprohackathon_G1 ubuntu@11.22.33.44:~

Then connect to the VM and update the system packages:
sudo apt update -y

Snakemake can be installed:
sudo apt install -y snakemake

Apptainer requires several dependencies:
sudo apt install -y build-essential libseccomp-dev pkg-config squashfs-tools cryptsetup uidmap

Then download and install Apptainer:
wget https://github.com/apptainer/apptainer/releases/download/v1.4.4/apptainer_1.4.4_amd64.deb
sudo dpkg -i apptainer_1.4.4_amd64.deb
sudo apt --fix-broken install -y

For reference, the versions used during development were Snakemake 7.32.4 and Apptainer 1.4.4.


----- Run the Pipeline 

Once all dependencies are installed, you can move to the project directory and launch the first part of the workflow with:
snakemake --use-singularity -j 4 -s Step1_RNAseq
This step performs the full RNA-seq processing workflow. Its execution time is approximately 2 hours on a 4-CPU IFB VM.

After it completes, you can run the second part of the analysis:
snakemake --use-singularity -j 4 -s Step2_Plots
This step runs the DESeq2 analysis, computes comparison statistics, and generates all figures. It is significantly faster, and results are produced in under 5 minutes.
