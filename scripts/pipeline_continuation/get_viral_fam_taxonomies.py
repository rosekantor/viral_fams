#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="Get taxonomies from ICTV genomes")

    parser.add_argument("-g", "--genomes-dir", help="genome sequences dir")
    parser.add_argument("-o", "--outfile", help="output")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit (-1)

    if args.genomes_dir is None:
        print ("must specify --{}\n".format('genomes-dir'))
        args_pass = False
    elif not os.path.exists(args.genomes_dir) or \
         not os.path.isdir(args.genomes_dir):
        print ("--{} {} must exist and not be empty\n".format('genomes-dir', args.genomes_dir))
        args_pass = False

    if args.outfile is None:
        print ("must specify --{}\n".format('outfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


def get_prot_tax(genomes_dir):
    tax_map = dict()

    for genome_root, genome_dirs, genome_files in os.walk(genomes_dir):
        for genome in genome_dirs:
                for info_root, info_dirs, info_files in os.walk(os.path.join(genome_root, genome)):
                    prot_ids = list()
                    taxonomy = ""
                    for file in info_files:
                        filepath = os.path.join(info_root, file)

                        if file.endswith("protein.faa"):
                            with open(filepath, "r") as file:
                                lines = file.readlines()
                                for line in lines:
                                    line = line.strip("*\n")
                                    if line.startswith(">"):
                                        prot_id = line.split(" ")[0][1:]
                                        prot_ids.append(prot_id)

                        elif file.endswith("genomic.gbff"):
                            for record in SeqIO.parse(filepath, "genbank"):
                                taxonomy = "; ".join(record.annotations["taxonomy"]) + f"; {record.annotations["organism"]}"    # using "_" instead of " " because hmmcan only uses first word
                    
                    if prot_ids and taxonomy:
                        for prot_id in prot_ids:
                            tax_map[prot_id] = taxonomy
                    else:
                        print(f"Did not find one of the files for {genome}\nprot_ids:{prot_ids}\ntaxonomy:{taxonomy}\n")

    return tax_map


def write_tax_to_file(tax_out, tax_map):
    with open(tax_out, "w") as file:
        for prot, tax in tax_map.items():
            file.write(f"{prot}\t{tax}\n")


def main():
    args = getargs()

    tax_map = get_prot_tax(args.genomes_dir)

    write_tax_to_file(args.outfile, tax_map)

    return 0


if __name__ == "__main__":
    main()