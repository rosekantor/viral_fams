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
    parser = argparse.ArgumentParser(description="Create test and training tables to be used for ml model.")

    parser.add_argument("-i", "--domain-table", help="domtblout output from hmmscan/hmmsearch")
    parser.add_argument("-m", "--hmmer-method", default="hmmscan", help="Options: hmmscan, hmmsearch")
    parser.add_argument("-c", "--tc-map", help="trusted cutoff map from determine_tc.py")
    parser.add_argument("-t", "--taxonomy-table", help="table that contains taxonomy of each gene from download_ICTV.ipynb")
    parser.add_argument("-n", "--nrows", type=int, default=None, required=False, help="Number of rows to randomly sampled. If none, all rows will be included.")
    parser.add_argument("--test-size", type=float, default=0, required=False, help="Proportion of rows that will be in the test dataset")
    parser.add_argument("-o", "--outdir", help="Output directory")
    
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
    
    if args.taxonomy_table is None:
        print ("must specify --{}\n".format('taxonomy-table'))
        args_pass = False
    elif not os.path.exists(args.taxonomy_table) or \
         not os.path.isfile(args.taxonomy_table) or \
         not os.path.getsize(args.taxonomy_table) > 0:
        print ("--{} {} must exist and not be empty\n".format('taxonomy-table', args.taxonomy_table))
        args_pass = False
    
    if args.test_size <= 0 or args.test_size >= 1:
        print ("--{} {} must be a between 0 and 1\n".format('test-size', args.test_size))
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


def read_taxonomy_table(filepath):
    df = pd.read_csv(filepath, sep="\t", usecols=["gene_id", "lineage"]).dropna()
    return dict(zip(df["gene_id"], df["lineage"]))


def identify_rank(rank_classification):
    if any([char in rank_classification for char in [" ", "-"]]):
        return "Species"
    for rank, suffixes in RANK_SUFFIX_MAP.items():
        for suffix in suffixes:
            if rank_classification.endswith(suffix):
                return rank
    return None


def create_encoding_map(hit_map, tax_map):
    all_rank_lineages = defaultdict(lambda: set([""]))
    all_rank_encodings = defaultdict(dict)
    gene_lineage_map = defaultdict(dict)

    all_genes = list(hit_map.keys())

    for gene in all_genes:
        lineage = tax_map[gene]
        split_lineage = lineage.split("; ")
        for i in range(len(split_lineage)):
            rank = identify_rank(split_lineage[i])
            if rank:
                rank_lineage = "; ".join(split_lineage[:i + 1])
                all_rank_lineages[rank].add(rank_lineage)

                gene_lineage_map[gene][rank] = rank_lineage

            else:
                print(f"Could not identify rank of '{split_lineage[i]}' in {lineage} from {gene}")
    
    for rank, lineages in all_rank_lineages.items():
        lineages = sorted(list(lineages))
        le = LabelEncoder()
        encodings = le.fit_transform(lineages)
        all_rank_encodings[rank] = dict(zip(lineages, encodings))
    
    return all_rank_encodings, gene_lineage_map


def write_encodings(outfile, all_rank_encodings):
    with open(outfile, "w") as file:
        file.write(f"Rank\tLineage\tEncoding\n")

        for rank in ORDERED_TAXA_RANKS:
            encoding_map = all_rank_encodings[rank]
            for lineage, encoding in encoding_map.items():
                file.write(f"{rank}\t{lineage}\t{encoding}\n")


def write_table(train_outfile, test_outfile, hit_map, tc_map, all_rank_encodings, gene_lineage_map, nrows, test_size):
    all_fams = sorted(list(tc_map.keys()))

    columns = ["gene"] + all_fams + ORDERED_TAXA_RANKS

    all_gene_rows = random.sample(list(hit_map.keys()), nrows) if nrows != None else list(hit_map.keys())

    n_test_rows = round(len(all_gene_rows) * test_size)
    test_rows = random.sample(all_gene_rows, n_test_rows)
    test_rows_set = set(test_rows)
    train_rows = [gene for gene in all_gene_rows if gene not in test_rows_set]

    for outfile, gene_rows in [(train_outfile, train_rows), (test_outfile, test_rows)]:
        with open(outfile, "w") as file:
            file.write(f"{'\t'.join(columns)}\n")

            for gene in gene_rows:
                hits = hit_map[gene]
                
                hit_info = [0] * len(all_fams)
                encoded_tax_info = [None] * len(ORDERED_TAXA_RANKS)

                for fam, bitscore, in hits.items():
                    hit_info[all_fams.index(fam)] = bitscore

                for i, rank in enumerate(ORDERED_TAXA_RANKS):
                    lineage = ""
                    if rank in gene_lineage_map[gene].keys():
                        lineage = gene_lineage_map[gene][rank]

                    encoded_tax_info[i] = all_rank_encodings[rank][lineage]
                
                file.write(f"{'\t'.join([str(x) for x in ([gene] + hit_info + encoded_tax_info)])}\n")

# main()
#
def main() -> int:
    args = getargs()

    hit_map = read_domtble_to_df(args.hmmer_method, args.domain_table)

    tc_map = read_tc_map(args.tc_map)

    tax_map = read_taxonomy_table(args.taxonomy_table)

    all_rank_encodings, gene_lineage_map = create_encoding_map(hit_map, tax_map)

    write_encodings(os.path.join(args.outdir, "encodings.tsv"), all_rank_encodings)

    write_table(os.path.join(args.outdir, "train.tsv"), os.path.join(args.outdir, "test.tsv"), hit_map, tc_map, all_rank_encodings, gene_lineage_map, args.nrows, args.test_size)
    
    print("Done")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())