#!/usr/bin/python3

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
    parser = argparse.ArgumentParser(description="Combine pulled genbank sequences into single files.")

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
    prot_seq_map = defaultdict(dict)
    nc_seq_map = defaultdict(dict)
    tax_map = dict()
    gff_map = defaultdict(list)

    for genome_root, genome_dirs, genome_files in os.walk(genomes_dir):
        for genome in genome_dirs:
                for info_root, info_dirs, info_files in os.walk(os.path.join(genome_root, genome)):
                    for file in info_files:
                        filepath = os.path.join(info_root, file)
                        
                        if file.endswith("protein.faa"):
                            fasta_dict = read_fasta_to_dict(filepath)
                            for prot, seq in fasta_dict.items():
                                prot_seq_map[genome][prot] = seq
                        
                        if file.endswith("genomic.fna") and not file.endswith("cds_from_genomic.fna"):
                            fasta_dict = read_fasta_to_dict(filepath)
                            for prot, seq in fasta_dict.items():
                                nc_seq_map[genome][prot] = seq

                        elif file.endswith("genomic.gbff"):
                            for record in SeqIO.parse(filepath, "genbank"):
                                taxonomy = "; ".join(record.annotations["taxonomy"]) + f"; {record.annotations["organism"]}"    # using "_" instead of " " because hmmcan only uses first word
                                tax_map[genome] = taxonomy

                        elif file.endswith("genomic.gff"):
                            with open(filepath) as gff:
                                lines = gff.readlines()
                                for line in lines:
                                    line = line.strip()
                                    if not line.startswith("#"):
                                        info = line.split("\t")
                                        if info[2] == "CDS":
                                            gff_map[genome].append(info)

    return prot_seq_map, nc_seq_map, tax_map, gff_map

# write_sequences_to_file(faa_out, seq_map, tax_map)
#
def write_sequences_to_file(faa_out, seq_map):
    with open(faa_out, "w") as file:
        for genome, prot_map in seq_map.items():
            for prot, seq in prot_map.items():
                file.write(f">{genome}~{prot}\n")
                file.write(f"{seq}\n")

def write_tax_to_file(tax_out, tax_map):
    with open(tax_out, "w") as file:
        for genome, tax in tax_map.items():
            file.write(f"{genome}\t{tax}\n")


def write_combined_gff(outdir, gff_map):
    gff_out = os.path.join(outdir, "sequences.gff")
    with open(gff_out, "w") as file:
        file.write("##gff-version  3\n")
        for genome, gff_infos in gff_map.items():
            for info in gff_infos:
                info[0] = f"{genome}"
                file.write(f"{'\t'.join(info)}\n")

# main()
#
def main() -> int:
    args = getargs()

    # get sequences and genome taxonomy
    # fam -> genome -> prot -> seq, genome -> tax, genome -> gff_info
    prot_seq_map, nc_seq_map, tax_map, gff_map = get_genome_info(args.genomes_dir)

    # write sequences
    faa_out = os.path.join(args.outdir, "sequences.faa")
    write_sequences_to_file(faa_out, prot_seq_map)

    fna_out = os.path.join(args.outdir, "sequences.fna")
    write_sequences_to_file(fna_out, nc_seq_map)

    tax_out = os.path.join(args.outdir, "taxonomy_map.tsv")
    write_tax_to_file(tax_out, tax_map)

    # write combined gff
    write_combined_gff(args.outdir, gff_map)

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