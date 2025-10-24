# Pauline : fasterq-dump via Homebrew
# Téléchargement des données en parallèle pour réduire le temps 

import os
import subprocess
import glob
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

SRA_IDS = [
    "SRR10379721", "SRR10379722", "SRR10379723",  # Persisters
    "SRR10379724", "SRR10379725", "SRR10379726"   # Controls
]
N_CPU = 4
N_PARALLEL = 3  # Nombre de téléchargements simultanés

# Création du dossier data
os.makedirs("data", exist_ok=True)
os.chdir("data")

def download_sra(SRA_ID):
    fq_gz = f"{SRA_ID}.fastq.gz"
    fq = f"{SRA_ID}.fastq"

    # Si fichier est déjà compressé on skip le téléchargement et la compression
    if os.path.exists(fq_gz):
        print(f"⏭️  Skipping {SRA_ID}: already downloaded & compressed.")
        return SRA_ID

    # Si fichier brut déjà présent , on skip le téléchargement 
    if os.path.exists(fq):
        print(f"⏭️  Skipping download for {SRA_ID}: .fastq already exists, will compress later.")
        return SRA_ID

    print(f"➡️  Downloading {SRA_ID} ...")
    try:
        subprocess.run([
            "fasterq-dump", "--threads", str(N_CPU), "--progress", SRA_ID
        ], check=True)
        print(f"✅ {SRA_ID} downloaded successfully.")
        return SRA_ID
    except subprocess.CalledProcessError:
        print(f"❌ Error downloading {SRA_ID}")
        return None


# Téléchargements parallèles
print(f"Starting parallel downloads ({N_PARALLEL} at a time)...")
with ThreadPoolExecutor(max_workers=N_PARALLEL) as executor:
    futures = {executor.submit(download_sra, sra): sra for sra in SRA_IDS}
    for future in as_completed(futures):
        sra = futures[future]
        if future.result() is None:
            print(f"⚠️ {sra} failed!")

# Compression (seulement pour les fichiers pas encore compressés)
fastq_files = glob.glob("*.fastq")
if fastq_files:
    print(" ℹ️ Compressing new FASTQ files...")
    subprocess.run("gzip *.fastq", shell=True, check=True)
else:
    print("ℹ️  No new FASTQ files to compress.")

# Vérification finale
gz_files = glob.glob("*.fastq.gz")
print("\nVerification:")
if len(gz_files) == len(SRA_IDS):
    print(f"✅ All {len(gz_files)} files are ready and compressed!")
else:
    missing = len(SRA_IDS) - len(gz_files)
    print(f"⚠️ Warning: {missing} file(s) missing!")
    print("Files found:", gz_files)
    sys.exit(1)
