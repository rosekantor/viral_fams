# download genome files from genbank ftp server

import sys
import os
from ftplib import FTP
import random
import tempfile
import gzip
import shutil
import argparse

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="Download a subset of species from genbank that are not in RefSeq")

    parser.add_argument("-n", "--num-species", type=int, help="Number of species to pull")
    parser.add_argument("-o", "--outdir", help="Output directory")
    
    args = parser.parse_args()
    args_pass = True

    if args.num_species is None:
        print ("must specify --{}\n".format('num-species'))
        args_pass = False

    if args.outdir is None:
        print ("must specify --{}\n".format('outdir'))
        args_pass = False
    elif os.path.exists(args.outdir):
        print(f"{args.outdir} already exists, please move or remove before running this script again.")
        args_pass = False
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args

def main():
    args = getargs()

    os.makedirs(args.outdir)

    print(f"Created directory: {args.outdir}")

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

    for i_species, species in enumerate(uncurated_species[0:args.num_species]):
        print(f"{i_species}. {species}")
        
        complete = False
        retries = 0
        while not complete:
            try:
                ftp.cwd(f'/genomes/genbank/viral/{species}/all_assembly_versions')
                for name, attributes in ftp.mlsd():
                    if attributes.get('type') == 'OS.unix=symlink':
                        genome = name
                        genome_outdir = os.path.join(args.outdir, genome)
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
                                        # print(f"Downloading {file}...")
                                        ftp.retrbinary(f'RETR {file}', f.write)

                                    unzipped_path = os.path.join(genome_outdir, file[:-len(".gz")])
                                    with gzip.open(temp_file, 'rb') as f_in, open(unzipped_path, 'wb') as f_out:
                                        # print(f"Unzipping {file}...")
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


if __name__ == '__main__':
    sys.exit(main())