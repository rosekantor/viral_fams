#!/usr/bin/python3
'''
Copyright 2024 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/gpl-3.0.txt>.
'''

import sys
import os
import argparse
import re
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
from Bio import SeqIO

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-g", "--genomes-dir", help="genome sequences dir")
    parser.add_argument("-o", "--outdir", help="output")
    
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

    if args.outdir is None:
        print ("must specify --{}\n".format('outdir'))
        args_pass = False
    elif not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    elif not os.path.isdir(args.outdir):
        print ("--{} {} is not a dir\n".format('outdir', args.outdir))

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args

# def get_genome_info(genomes_dir)
#
def get_genome_info(genomes_dir):
    seq_map = defaultdict(dict)
    tax_map = dict()
    for genome_root, genome_dirs, genome_files in os.walk(genomes_dir):
        for genome in genome_dirs:
                for info_root, info_dirs, info_files in os.walk(os.path.join(genome_root, genome)):
                    for file in info_files:
                        filepath = os.path.join(info_root, file)
                        
                        if file.endswith("protein.faa"):
                            fasta_dict = read_fasta_to_dict(filepath)
                            for prot, seq in fasta_dict.items():
                                seq_map[genome][prot] = seq

                        elif file.endswith("genomic.gbff"):
                            for record in SeqIO.parse(filepath, "genbank"):
                                taxonomy = "; ".join(record.annotations["taxonomy"]) + f"; {record.annotations["organism"]}"    # using "_" instead of " " because hmmcan only uses first word
                                tax_map[genome] = taxonomy

    return seq_map, tax_map

# write_fam_sequences_to_faa(outdir, seq_map, tax_map)
#
def write_fam_sequences_to_faa(outdir, seq_map, tax_map):
    faa_out = os.path.join(outdir, "sequences.faa")
    with open(faa_out, "w") as file:
        for genome, prot_map in seq_map.items():
            for prot, seq in prot_map.items():
                file.write(f">{genome}_{prot}\n")
                file.write(f"{seq}\n")
    
    tax_out = os.path.join(outdir, "taxonomy_map.tsv")
    with open(tax_out, "w") as file:
        for genome, tax in tax_map.items():
            file.write(f"{genome}\t{tax}\n")

# main()
#
def main() -> int:
    args = getargs()

    # get sequences and genome taxonomy
    # fam -> genome -> prot -> seq, genome -> tax
    seq_map, tax_map = get_genome_info(args.genomes_dir)

    # write fam sequences
    write_fam_sequences_to_faa(args.outdir, seq_map, tax_map)

    print ("DONE")
    return 0

# read_fasta_to_dict(filepath)
#
def read_fasta_to_dict(filepath):
    seq_map = dict()
    with open(filepath, "r") as file:
        lines = file.readlines()
        name = ""
        seq = ""
        for line in lines:
            line = line.strip("*\n")
            if line.startswith(">"):
                if name and seq:
                    seq_map[name] = seq
                name = line.split(" ")[0][1:]
                seq = ""
            else:
                seq += line
                
        if name and seq:
            seq_map[name] = seq
    
    return seq_map

# exec()
#
if __name__ == '__main__':
    sys.exit(main())