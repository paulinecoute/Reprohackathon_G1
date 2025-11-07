# Snakefile 
# mis à jour le 06.11
# les images docker doivent avoir été construites en amont


SRR = ["SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726"]

# Dossiers pour mettre les données au fur et à mesure de l'analyse
raw_directory = "data_raw"
trimmed_directory = "data_trimmed"
reference_directory = "reference"
mapping_directory = "data_mapped"
count_directory = "data_counts"

rule all:
    input:
        f"{count_directory}/counts.txt"

# Première étape : Téléchargement des fastq
# Utilisation de Dockerfile

rule download_and_compress:
    output: f"{raw_directory}/{{srr}}.fastq.gz"
    params: threads = 4
    wildcard_constraints: srr="SRR[0-9]+"
    container: "docker:base"
    shell:
        r"""
        mkdir -p {raw_directory}
        cd {raw_directory}
        fasterq-dump --threads {params.threads} --progress {wildcards.srr}
        pigz -p {params.threads} {wildcards.srr}.fastq
        """

# Deuxième étape : Trimming
# Utilisation de Dockerfile.trim        

rule trim_reads:
    input:
        fastq = f"{raw_directory}/{{srr}}.fastq.gz"
    output:
        f"{trimmed_directory}/{{srr}}_trimmed.fq.gz"
    params:
        quality = 20, length = 25
    wildcard_constraints: srr="SRR[0-9]+"
    container: "docker:trim"
    shell:
        r"""
        mkdir -p {trimmed_directory}
        cd {trimmed_directory}
        trim_galore -q {params.quality} --phred33 --length {params.length} ../{input.fastq}
        """

# Troisème étape : téléchargement du génome de référence et de ses annotations
# De base je ne voulais pas utiliser de dockerfile mais sur mon pc wget ne fonctionne pas 

rule download_reference:
    output:
        fasta = f"{reference_directory}/reference.fasta",
        gff = f"{reference_directory}/reference.gff"
    container: "docker:base"
    shell:
        r"""
        mkdir -p {reference_directory}
        cd {reference_directory}
        wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta&retmode=text"
        wget -q -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=gff3&format=text"
        """

# Quatrième étape : indexaction
# Docker.bowtie

rule build_index:
    input:
        fasta = f"{reference_directory}/reference.fasta"
    output:
        f"{reference_directory}/reference.1.bt2"
    container: "docker:bowtie"
    shell:
        "bowtie-build {input.fasta} {reference_directory}/reference"

# Cinquième étape : mapping
# Docker.bowtie (c'est le même que pour celui d'avant)
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
    wildcard_constraints: srr="SRR[0-9]+"
    container: "docker:bowtie"
    shell:
        r"""
        mkdir -p {mapping_directory}
        bowtie -p {params.threads} -S {params.index_name} <(gunzip -c {input.fastq}) | \
            samtools sort -@ {params.threads} -o {output.bam}
        samtools index {output.bam}
        """

# Sixième étape : comptage des reads 
# Docker.featurecounts

rule count_reads:
    input:
        gff = f"{reference_directory}/reference.gff",
        bam = expand(f"{mapping_directory}/{{srr}}.bam", srr=SRR)
    output:
        f"{count_directory}/counts.txt"
    params: threads = 4
    container: "docker:featurecounts"
    shell:
        r"""
        mkdir -p {count_directory}
        featureCounts -t gene -g ID -F GTF -T {params.threads} \
            -a {input.gff} -o {output} {mapping_directory}/*.bam
        """
