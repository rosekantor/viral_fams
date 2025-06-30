# download genome files from genbank ftp server

import sys
import os
from ftplib import FTP
import random
import tempfile
import gzip
import shutil

outdir = "/p/vast1/ibap/chivian1/proj/viral_fams/data/uncurated_genbank/genomes"

try:
    ftp = FTP('ftp.ncbi.nlm.nih.gov', user="anonymous", passwd="golez1@llnl.gov")
    print("Successfully connected and logged in to FTP server.")
except Exception as e:
    print(f"Error connecting to FTP server: {e}")
    sys.exit(-1)

ftp.cwd('/genomes/refseq/viral')

refseq_species = ftp.nlst()

ftp.cwd('/genomes/genbank/viral')

genbank_species = ftp.nlst()

uncurated_species = [species for species in genbank_species if species not in refseq_species]

# uncurated_species.sort()
random.shuffle(uncurated_species)

n_species = 50
for species in uncurated_species[0:50]:
    print(species)
    
    complete = False
    retries = 0
    while not complete:
        try:
            ftp.cwd(f'/genomes/genbank/viral/{species}/all_assembly_versions')
            for genome in ftp.nlst():
                genome_outdir = os.path.join(outdir, genome)
                os.makedirs(genome_outdir, exist_ok=True)

                ftp.cwd(f'/genomes/genbank/viral/{species}/all_assembly_versions/{genome}')
                for file in ftp.nlst():
                    if any(file.endswith(suffix) for suffix in ["genomic.fna.gz", "genomic.gbff.gz", "genomic.gff.gz", "genomic.gtf.gtf", "protein.faa.gz"]):
                        unzipped_path = os.path.join(genome_outdir, file[:-len(".gz")])

                        if os.path.exists(unzipped_path) and \
                        os.path.isfile(unzipped_path) and \
                        os.path.getsize(unzipped_path) > 0:
                            continue

                        with tempfile.TemporaryDirectory() as tmpdir:
                            temp_file = os.path.join(tmpdir, file)
                            with open(temp_file, 'wb') as f:
                                print(f"Downloading {file}...")
                                ftp.retrbinary(f'RETR {file}', f.write)

                            unzipped_path = os.path.join(genome_outdir, file[:-len(".gz")])
                            with gzip.open(temp_file, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
                                print(f"Unzipping {file}...")
                                shutil.copyfileobj(f_in, f_out)
            
            complete = True

        except Exception as e1:
            retries += 1
            if retries > 10:
                print(f"{e1}\nToo many retries... Exiting.")
                sys.exit(-1)
            else:
                print(f"{e1}\nTrying again...")

                try:
                    ftp = FTP('ftp.ncbi.nlm.nih.gov', user="anonymous", passwd="golez1@llnl.gov")
                    print("Successfully connected and logged in to FTP server.")
                except Exception as e2:
                    print(f"Error connecting to FTP server: {e2}")
                    sys.exit(-1)


