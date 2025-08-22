#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
import random
from sklearn.preprocessing import LabelEncoder

RANK_SUFFIX_MAP = {
    "Domain": ["Viruses"],
    "Realm": ["viria"],
    "Kingdom": ["virae"],
    "Phylum": ["viricota"],
    "Subphylum": ["viricotina"],
    "Class": ["viricetes"],
    "Order": ["virales"],
    "Suborder": ["virineae"],
    "Family": ["viridae", "formidae", "satellitidae"],
    "Subfamily": ["virinae", "satellitinae"],
    "Genus": ["virus", "form", "satellite"]
}
# note: Species can be identified by having a space

ORDERED_TAXA_RANKS = ["Domain", "Realm", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Suborder", "Family", "Subfamily", "Genus", "Species"]

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="Create novel table that will be inputted into ML model.")

    parser.add_argument("-i", "--domain-table", help="domtblout output from hmmscan/hmmsearch")
    parser.add_argument("-m", "--hmmer-method", default="hmmscan", help="Options: hmmscan, hmmsearch")
    parser.add_argument("-c", "--tc-map", help="trusted cutoff map from determine_tc.py")
    parser.add_argument("-o", "--outdir", help="output directory")
    
    args = parser.parse_args()
    args_pass = True

    if args.domain_table is None:
        print ("must specify --{}\n".format('domain-table'))
        args_pass = False
    if not os.path.isfile(args.domain_table) or \
         not os.path.isfile(args.domain_table) or \
         not os.path.getsize(args.domain_table) > 0:
        print ("--{} {} must exist and not be empty\n".format('domain-table', args.domain_table))
        args_pass = False

    if args.hmmer_method not in ["hmmscan", "hmmsearch"]:
        print ("--{} must be a valid option\n".format('hmmer-method'))
        args_pass = False
    
    if args.tc_map is None:
        print ("must specify --{}\n".format('tc-map'))
        args_pass = False
    elif not os.path.exists(args.tc_map) or \
         not os.path.isfile(args.tc_map) or \
         not os.path.getsize(args.tc_map) > 0:
        print ("--{} {} must exist and not be empty\n".format('tc-map', args.tc_map))
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


# read_domtble_to_df()
#
def read_domtble_to_df(hmmer_method, filepath):
    hit_map = defaultdict(dict)

    if hmmer_method == "hmmscan":
        type = "hmmscan3-domtab"
    elif hmmer_method == "hmmsearch":
        type = "hmmsearch3-domtab"

    with open(filepath, 'r') as file:
        for qresult in SearchIO.parse(file, type):
            for hit in qresult.hits:
                if hmmer_method == "hmmscan":
                    gene = qresult.id
                    fam = hit.id
                elif hmmer_method == "hmmsearch":
                    gene = hit.id
                    fam = qresult.id
                bitscore = hit.bitscore

                hit_map[gene][fam] = bitscore

    return hit_map


# read_tc_map()
#
def read_tc_map(filepath):
    map = dict()
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            map[info[0]] = float(info[1])
    
    return map


def write_table(outfile, hit_map, tc_map):
    all_fams = sorted(list(tc_map.keys()))

    columns = ["gene"] + all_fams

    with open(outfile, "w") as file:
        file.write(f"{'\t'.join(columns)}\n")

        for gene in list(hit_map.keys()):
            hits = hit_map[gene]
            
            hit_info = [0] * len(all_fams)

            for fam, bitscore, in hits.items():
                hit_info[all_fams.index(fam)] = bitscore
            
            file.write(f"{'\t'.join([str(x) for x in ([gene] + hit_info)])}\n")

# main()
#
def main() -> int:
    args = getargs()

    hit_map = read_domtble_to_df(args.hmmer_method, args.domain_table)

    tc_map = read_tc_map(args.tc_map)

    write_table(os.path.join(args.outdir, "input_table.tsv"), hit_map, tc_map)
    
    print("Done")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())