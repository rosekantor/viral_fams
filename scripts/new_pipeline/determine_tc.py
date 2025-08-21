#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--domain-table", help="")
    parser.add_argument("-m", "--hmmer-method", default="hmmscan", help="Options: hmmscan, hmmsearch")
    parser.add_argument("-f", "--fam-members-file", help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.domain_table is None:
        print ("must specify --{}\n".format('domain-table'))
        args_pass = False
    if not os.path.isfile(args.domain_table):
        print ("--{} {} must exist and not be empty\n".format('domain-table', args.domain_table))
        args_pass = False

    if args.hmmer_method not in ["hmmscan", "hmmsearch"]:
        print ("--{} must be a valid option\n".format('hmmer-method'))
        args_pass = False

    if args.fam_members_file is None:
        print ("must specify --{}\n".format('fam-members-file'))
        args_pass = False
    elif not os.path.exists(args.fam_members_file) or \
         not os.path.isfile(args.fam_members_file) or \
         not os.path.getsize(args.fam_members_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('fam-members-file', args.fam_members_file))
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
    hit_bitscores = defaultdict(lambda: defaultdict(list))

    with open(filepath, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                fields = line.split()[0:22]
                if hmmer_method == "hmmscan":
                    fam = fields[0]
                    gene = fields[3]
                elif hmmer_method == "hmmsearch":
                    gene = fields[0]
                    fam = fields[3]

                # score = fields[13]
                score = fields[7]

                hit_bitscores[fam][gene].append(score)

    return hit_bitscores


def read_fam_members(filepath):
    members_map = defaultdict(list)

    df = pd.read_csv(filepath, sep="\t", usecols=["#GroupName", "ProteinIDs"])

    for i, row in df.iterrows():
        group = row["#GroupName"]
        protein_ids = row["ProteinIDs"].split(",")

        members_map[group] = protein_ids

    return members_map


def determine_fam_tcs(hit_bitscores, fam_prot_ids):
    fam_tc_map = dict()
    
    for fam, prot_ids in fam_prot_ids.items():
        all_bitscores = [0]
        for prot in prot_ids:
            if prot in hit_bitscores[fam].keys():
                all_bitscores += hit_bitscores[fam][prot]
        
        fam_tc_map[fam] = min(all_bitscores)
    
    return fam_tc_map


def write_dict_to_file(outfile, map):
    with open(outfile, "w") as file:
        for key, val in map.items():
            file.write(f"{key}\t{val}\n")


def main():
    args = getargs()

    hit_bitscores = read_domtble_to_df(args.hmmer_method, args.domain_table)

    fam_prot_ids = read_fam_members(args.fam_members_file)

    fam_tc_map = determine_fam_tcs(hit_bitscores, fam_prot_ids)

    tc_outfile = os.path.join(args.outdir, "tc_map.tsv")
    write_dict_to_file(tc_outfile, fam_tc_map)

    return 0


if __name__ == "__main__":
    main()