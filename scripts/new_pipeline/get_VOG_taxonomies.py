#!/usr/bin/python3

# This script isn't needed because VOGDB has a lca file


import sys
import os
import argparse
import pandas as pd
from collections import defaultdict

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-m", "--fam-members-file", help="")
    parser.add_argument("-t", "--prot-tax-map", help="")
    parser.add_argument("-o", "--outfile", help="output")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit (-1)

    if args.fam_members_file is None:
        print ("must specify --{}\n".format('fam-members-file'))
        args_pass = False
    elif not os.path.exists(args.fam_members_file) or \
         not os.path.isfile(args.fam_members_file) or \
         not os.path.getsize(args.fam_members_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('fam-members-file', args.fam_members_file))
        args_pass = False
    
    if args.prot_tax_map is None:
        print ("must specify --{}\n".format('prot-tax-map'))
        args_pass = False
    elif not os.path.exists(args.prot_tax_map) or \
         not os.path.isfile(args.prot_tax_map) or \
         not os.path.getsize(args.prot_tax_map) > 0:
        print ("--{} {} must exist and not be empty\n".format('prot-tax-map', args.prot_tax_map))
        args_pass = False

    if args.outfile is None:
        print ("must specify --{}\n".format('outfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


def read_fam_members(filepath):
    members_map = defaultdict(list)

    df = pd.read_csv(filepath, sep="\t", usecols=["#GroupName", "ProteinIDs"])

    for i, row in df.iterrows():
        group = row["#GroupName"]
        protein_ids = row["ProteinIDs"].split(",")

        for j, protein in enumerate(protein_ids):
            protein_ids[j] = ".".join(protein.split(".")[1:])

        members_map[group] = protein_ids

    return members_map


def read_tax_map(filepath):
    tax_map = dict()

    with open(filepath) as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            tax_map[info[0]] = info[1]

    return tax_map


def find_lca_taxa(all_taxa):
    all_split_taxa = [taxa.split("; ") for taxa in all_taxa]
    lca_index = 0
    if len(all_split_taxa) > 0:
        while True:
            if not all(lca_index < len(taxa) for taxa in all_split_taxa):
                break
            if not all(taxa[lca_index] == all_split_taxa[0][lca_index] for taxa in all_split_taxa):
                break

            lca_index += 1

        return "; ".join(all_split_taxa[0][:lca_index])
    
    return ""


def get_vog_tax_map(vog_prot_ids, prot_tax_map):
    vog_tax_map = dict()

    for vog, prot_ids in vog_prot_ids.items():
        taxonomy = find_lca_taxa([prot_tax_map[prot_id] for prot_id in prot_ids if prot_id in prot_tax_map.keys()])
        vog_tax_map[vog] = taxonomy

    return vog_tax_map


def write_tax_to_file(tax_out, tax_map):
    with open(tax_out, "w") as file:
        for prot, tax in tax_map.items():
            file.write(f"{prot}\t{tax}\n")


def main():
    args = getargs()

    vog_prot_ids = read_fam_members(args.fam_members_file)

    prot_tax_map = read_tax_map(args.prot_tax_map)

    vog_tax_map = get_vog_tax_map(vog_prot_ids, prot_tax_map)

    write_tax_to_file(args.outfile, vog_tax_map)

    return 0


if __name__ == "__main__":
    main()