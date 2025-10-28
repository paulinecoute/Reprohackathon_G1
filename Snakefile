# Étape 0 

# Fichiers FASTQ (run)
SRR = ["SRR10379721", "SRR10379722", "SRR10379723","SRR10379724", "SRR10379725", "SRR10379726"]

# Dossiers de travail 
raw_directory = "data_raw"  
trimmed_directory = "data_trimmed"
reference_directory = "reference"
mapping_directory = "data_mapped"
count_directory = "data_counts"

rule all:
    input:
        f"{count_directory}/counts.txt"

# Étape 1 : Téléchargement et compression des fichiers fastq 
# dockerfile sra-toolkits

rule download_and_compress:
    output:
        f"{raw_directory}/{{srr}}.fastq.gz"
    params:
        threads = 4
    shell:
        r"""
        mkdir -p {raw_directory}
        cd {raw_directory}
        fasterq-dump --threads {params.threads} --progress {wildcards.srr}
        pigz -p {params.threads} {wildcards.srr}.fastq
        echo "{wildcards.srr}.fastq.gz downloaded and compressed"
        """

# Étape 2 : Trimming
# dockerfile trimming

rule trim_reads:
    input:
        f"{raw_directory}/{{srr}}.fastq.gz"
    output:
        f"{trimmed_directory}/{{srr}}_trimmed.fq.gz"
    params:
        quality = 20,
        length = 25
    shell:
        r"""
        mkdir -p {trimmed_directory}
        cd {trimmed_directory}
        trim_galore -q {params.quality} --phred33 --length {params.length} ../{input}
        echo "{wildcards.srr} trimmed"
        """

# Étape 3 : Téléchargement du génome de référence et annotations
# Pas de dockerfile nécessaire, il faut simplement avoir wget sur son pc

rule download_reference:
    output:
        f"{reference_directory}/reference.fasta",
        f"{reference_directory}/reference.gff"
    shell:
        r"""
        mkdir -p {reference_directory}
        cd {reference_directory}
        wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        wget -q -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
        echo "Reference genome and annotations downloaded"
        """

# Étape 4 : Création de l'index du génome
# DOCKERFILE A CREER

rule build_index:
    input:
        f"{reference_directory}/reference.fasta"
    output:
        f"{reference_directory}/reference.1.bt2"
    params:
        index_name = f"{reference_directory}/reference"
    shell:
        r"""
        cd {reference_directory}
        bowtie-build reference.fasta reference
        echo "Genome index created"
        """

# Étape 5 : Mapping
# DOCKERFILE A CREER

rule mapping:
    input:
        fastq = f"{trimmed_directory}/{{srr}}_trimmed.fq.gz",
        index = f"{reference_directory}/reference.1.bt2"
    output:
        bam = f"{mapping_directory}/{{srr}}.bam",
        bai = f"{mapping_directory}/{{srr}}.bam.bai"
    params:
        threads = 4,
        index_name = f"{reference_directory}/reference"
    shell:
        r"""
        mkdir -p {mapping_directory}
        cd {mapping_directory}
        bowtie -p {params.threads} -S {params.index_name} <(gunzip -c ../{input.fastq}) | \
            samtools sort -@ {params.threads} -o {wildcards.srr}.bam
        samtools index {wildcards.srr}.bam
        echo "{wildcards.srr} mapped and indexed"
        """

# Étape 6 : Comptage des reads
# DOCKERFILE A CREER

rule count_reads:
    input:
        gff = f"{reference_directory}/reference.gff",
        bam = expand(f"{mapping_directory}/{{srr}}.bam", srr=SRR)
    output:
        f"{count_directory}/counts.txt"
    params:
        threads = 4
    shell:
        r"""
        mkdir -p {count_directory}
        cd {count_directory}
        featureCounts --extraAttributes Name -t gene -g ID -F GTF \
            -T {params.threads} -a ../{input.gff} \
            -o counts.txt ../{mapping_directory}/*.bam
        echo "Read counting finished, results saved in {count_directory}/counts.txt"
        """


