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

    parser.add_argument("-i", "--lin-hits-table", help="")
    parser.add_argument("-a", "--hmm-hits-dir", help="genome hits dir")
    parser.add_argument("-d", "--db-name", help="")
    parser.add_argument("-g", "--genomes-dir", help="genome sequences dir")
    parser.add_argument("-o", "--outdir", help="output dir")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 11:
        parser.print_help()
        sys.exit (-1)

    if args.lin_hits_table is None:
        print ("must specify --{}\n".format('lin-hits-table'))
        args_pass = False
    elif not os.path.exists(args.lin_hits_table) or \
         not os.path.isfile(args.lin_hits_table) or \
         not os.path.getsize(args.lin_hits_table) > 0:
        print ("--{} {} must exist and not be empty\n".format('lin-hits-table', args.lin_hits_table))
        args_pass = False

    if args.hmm_hits_dir is None:
        print ("must specify --{}\n".format('hmm-hits-dir'))
        args_pass = False
    elif not os.path.exists(args.hmm_hits_dir) or \
         not os.path.isdir(args.hmm_hits_dir):
        print ("--{} {} must exist and not be empty\n".format('hmm-hits-dir', args.hmm_hits_dir))
        args_pass = False

    if args.genomes_dir is None:
        print ("must specify --{}\n".format('genomes-dir'))
        args_pass = False
    elif not os.path.exists(args.genomes_dir) or \
         not os.path.isdir(args.genomes_dir):
        print ("--{} {} must exist and not be empty\n".format('genomes-dir', args.genomes_dir))
        args_pass = False

    if args.db_name is None:
        print ("must specify --{}\n".format('db-name'))
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

# read_lin_hits_to_dict(filepath)
#
def read_lin_hits_to_dict(filepath):
    prot_fam_map = defaultdict(lambda: defaultdict(list))
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    genome_hits_map = {fam: row["GENOME_HITS"] for fam, row in df.iterrows()}

    for fam, genomes in genome_hits_map.items():
        genomes = genomes.split(";")
        for genome in genomes:
            genome = genome.split(":")
            genome_id = genome[0]
            prot_id = genome[1]
            prot_fam_map[genome_id][prot_id].append(fam)

    return prot_fam_map

# get_hit_coords(hmm_hits_dir, db_name, prot_fam_map)
#
def get_hit_coords(hmm_hits_dir, db_name, prot_fam_map):
    prot_coords = defaultdict(lambda: defaultdict(dict))
    for annot_root, annot_dirs, annot_files in os.walk(hmm_hits_dir):
        for genome in annot_dirs:
            if genome in prot_fam_map.keys():
                for genome_root, genome_dirs, genome_files in os.walk(os.path.join(annot_root, genome)):
                    for file in genome_files:
                        if db_name in file and "domtbl" in file:
                            seen_prot_hits = set()
                            with open(os.path.join(genome_root, file), 'r') as file:
                                for qresult in SearchIO.parse(file, 'hmmscan3-domtab'):
                                    for hit in qresult.hits:
                                        for hsp in hit.hsps:
                                            query = qresult.id
                                            target = hit.id
                                            env_start = hsp.env_start
                                            env_end = hsp.env_end
                                            if (query, target) not in seen_prot_hits:
                                                if query in prot_fam_map[genome].keys():
                                                    if target in prot_fam_map[genome][query]:
                                                        prot_coords[genome][query][target] = (env_start, env_end)
                                                        seen_prot_hits.add((query, target))
    
    return prot_coords

# def get_genome_info(genomes_dir, prot_coords)
#
def get_genome_info(genomes_dir, prot_coords):
    seq_map = defaultdict(lambda: defaultdict(dict))
    tax_map = dict()
    for genome_root, genome_dirs, genome_files in os.walk(genomes_dir):
        for genome in genome_dirs:
            if genome in prot_coords.keys():
                for info_root, info_dirs, info_files in os.walk(os.path.join(genome_root, genome)):
                    for file in info_files:
                        filepath = os.path.join(info_root, file)
                        
                        if file.endswith("-protein.faa"):
                            fasta_dict = read_fasta_to_dict(filepath)
                            for prot, seq in fasta_dict.items():
                                if prot in prot_coords[genome].keys():
                                    for fam, coords in prot_coords[genome][prot].items():
                                        seq_map[fam][genome][prot] = seq[coords[0]:coords[1]]

                        elif file.endswith("-genomic.gbff"):
                            for record in SeqIO.parse(filepath, "genbank"):
                                taxonomy = "; ".join(record.annotations["taxonomy"]) + f"; {record.annotations["organism"]}"
                                tax_map[genome] = taxonomy

    return seq_map, tax_map

# write_fam_sequences_to_faa(outdir, seq_map, tax_map)
#
def write_fam_sequences_to_faa(outdir, seq_map, tax_map):
    for fam, genome_map in seq_map.items():
        with open(os.path.join(outdir, f"{fam}-sliced.faa"), "w") as file:
            for genome, prot_map in genome_map.items():
                for prot, seq in prot_map.items():
                    file.write(f">{genome}:{prot}:{tax_map[genome]}\n")
                    file.write(f"{seq}\n")

# main()
#
def main() -> int:
    args = getargs()

    # get fam prot hits
    # genome -> protein -> fams_list
    prot_fam_map = read_lin_hits_to_dict(args.lin_hits_table)

    # get hit coords
    # genome -> protein -> fam: (start, end)
    prot_coords = get_hit_coords(args.hmm_hits_dir, args.db_name, prot_fam_map)

    # get fam sequences and genome taxonomy
    # fam -> genome -> prot -> seq, genome -> tax
    seq_map, tax_map = get_genome_info(args.genomes_dir, prot_coords)

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

