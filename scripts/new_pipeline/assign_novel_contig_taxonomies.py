#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.pyplot as plt


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--domain-table", help="")
    parser.add_argument("-m", "--hmmer-method", default="hmmscan", help="Options: hmmscan, hmmsearch")
    parser.add_argument("-c", "--tc-map", help="")
    parser.add_argument("-t", "--fam-tax-map", help="")
    parser.add_argument("-s", "--fam-host-file", help="")
    parser.add_argument("-o", "--outdir", help="")
    
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

    if args.fam_tax_map is None:
        print ("must specify --{}\n".format('fam-tax-map'))
        args_pass = False
    elif not os.path.exists(args.fam_tax_map) or \
         not os.path.isfile(args.fam_tax_map) or \
         not os.path.getsize(args.fam_tax_map) > 0:
        print ("--{} {} must exist and not be empty\n".format('fam-tax-map', args.fam_tax_map))
        args_pass = False

    if args.fam_host_file is None:
        print ("must specify --{}\n".format('fam-host-file'))
        args_pass = False
    elif not os.path.exists(args.fam_host_file) or \
         not os.path.isfile(args.fam_host_file) or \
         not os.path.getsize(args.fam_host_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('fam-host-file', args.fam_host_file))
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
    hits = list()

    if hmmer_method == "hmmscan":
        type = "hmmscan3-domtab"
    elif hmmer_method == "hmmsearch":
        type = "hmmsearch3-domtab"

    with open(filepath, 'r') as file:
        for qresult in SearchIO.parse(file, type):
            for hit in qresult.hits:
                for hsp in hit.hsps:
                    if hmmer_method == "hmmscan":
                        gene = qresult.id
                        fam = hit.id
                    elif hmmer_method == "hmmsearch":
                        gene = hit.id
                        fam = qresult.id
                    bitscore = hsp.bitscore
                    # bitscore = hit.bitscore
                    evalue = hsp.evalue
                    env_start = hsp.env_start
                    env_end = hsp.env_end

                    hits.append((gene, fam, bitscore, evalue, env_start, env_end))

    return pd.DataFrame(hits, columns=["gene", "fam", "score", "evalue", "env start", "env end"])


# read_fasta_to_dict()
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
                name = line.split(" # ")[0][1:]
                seq = ""
            else:
                seq += line
                
        if name and seq:
            seq_map[name] = seq
    
    return seq_map

# read_tax_map()
#
def read_tax_map(filepath):
    df = pd.read_csv(filepath, sep="\t", usecols=["#GroupName", "LastCommonAncestor_Name"]).fillna("")

    df["LastCommonAncestor_Name"] = df["LastCommonAncestor_Name"].str.replace(";", "; ")

    return dict(zip(df["#GroupName"], df["LastCommonAncestor_Name"]))


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


# read_fam_host_file()
#
def read_fam_host_file(filepath):
    df = pd.read_csv(filepath, sep="\t", usecols=["Group", "Best Host"], keep_default_na=False)
    return dict(zip(df["Group"], df["Best Host"]))


# assign_gene_taxonomy()
#
def assign_gene_taxonomy(domtbl_df, tc_map):
    gene_fam_map = dict()
    gene_bitscore_map = dict()

    gene_grouped_tbl_df = domtbl_df.groupby("gene")
    for gene, df in gene_grouped_tbl_df:
        df.sort_values(by="score", ascending=False) # get largest bitscores first
        for i, info in df.iterrows():
            fam = info["fam"]
            bitscore = info["score"]
            # evalue = info["evalue"]
            # env_start = info["env start"]
            # env_end = info["env end"]

            if bitscore >= tc_map[fam]:
                gene_fam_map[gene] = fam
                gene_bitscore_map[gene] = bitscore
                break
    
    return gene_fam_map, gene_bitscore_map


# write_gene_tax_to_file()
#
def write_gene_tax_to_file(tax_out, tax_map, bitscore_map):
    with open(tax_out, "w") as file:
        for gene, tax in tax_map.items():
            file.write(f"{gene}\t{bitscore_map[gene]}\t{tax}\n")

# write_contig_tax_to_file()
#
def write_contig_tax_to_file(tax_out, tax_map, host_map):
    with open(tax_out, "w") as file:
        file.write("Genome\tTaxonomy\tHost\n")
        for genome, tax in tax_map.items():
            file.write(f"{genome}\t{tax}\t{host_map[genome]}\n")


# get_contig_id()
#
def get_contig_id(query_name):
    # return "_".join(query_name.split("_")[:-1])
    return query_name.split("~")[0]


# get_contigs()
#
def get_contigs(queries):
    contig_all_queries_map = defaultdict(list)
    for query in queries:
        contig_id = get_contig_id(query)
        contig_all_queries_map[contig_id].append(query)
    
    return contig_all_queries_map


# get_contig_host_counts()
#
def get_contig_host_counts(contig_prots_map, gene_fam_map, fam_host_map):
    contig_host_counts = defaultdict(lambda: defaultdict(int))
    for contig, prots in contig_prots_map.items():
        for prot in prots:
            fam = gene_fam_map[prot]
            best_fam_host = fam_host_map[fam]

            contig_host_counts[contig][best_fam_host] += 1
    
    contig_host_map = dict()
    for contig, host_counts in contig_host_counts.items():
        contig_host_map[contig] = max(host_counts, key=host_counts.get)

    return contig_host_counts, contig_host_map


# write_contig_host_counts()
#
def write_contig_host_counts(outfile, contig_host_counts, contig_host_map):
    all_hosts = sorted(list(set(host for host_counts in contig_host_counts.values() for host in host_counts.keys())))

    with open(outfile, "w") as file:
        file.write("Contig\tBest Host")
        for host in all_hosts:
            file.write(f"\t{host}")
        file.write("\n")

        for contig, host_counts in contig_host_counts.items():
            best_host = contig_host_map[contig]

            file.write(f"{contig}\t{best_host}")
            for host in all_hosts:
                if host in host_counts.keys():
                    file.write(f"\t{host_counts[host]}")
                else:
                    file.write(f"\t{0}")
            file.write("\n")


# get_greatest_keys()
#
def get_greatest_keys(counts_dict):
    largest_keys = list()
    largest_count = 0

    for key, count in counts_dict.items():
        if count > largest_count:
            largest_keys = [key]
            largest_count = count
        elif count == largest_count:
            largest_keys.append(key)

    return largest_keys


# find_lca_taxa()
#
def find_lca_taxa(all_split_taxa):
    lca_index = 0
    while True:
        if not all(lca_index < len(taxa) for taxa in all_split_taxa):
            break
        if not all(taxa[lca_index] == all_split_taxa[0][lca_index] for taxa in all_split_taxa):
            break
        lca_index += 1

    return "; ".join(all_split_taxa[0][:lca_index])


# decide_most_common_taxonomy()
#
def decide_most_common_taxonomy(all_split_taxa):
    level_idx = 0

    while True:
        # remove taxa not in level
        all_split_taxa = [taxa for taxa in all_split_taxa if level_idx < len(taxa)]

        if len(all_split_taxa) > 0:
            # get most common rank classification
            rank_counts = defaultdict(int)
            for taxa in all_split_taxa:
                rank_counts[taxa[level_idx]] += 1

            most_common_rank_classifications = get_greatest_keys(rank_counts)

            # remove taxa that do not have any of the most common classifications
            remove_indexes = list()
            for i, taxa in enumerate(all_split_taxa):
                if not any([taxa[level_idx] == classification  for classification in most_common_rank_classifications]):
                    remove_indexes.append(i)

            for i, j in enumerate(remove_indexes):
                remove_index = j - i
                all_split_taxa.pop(remove_index)
            
            # get lca if there is a tie
            if len(most_common_rank_classifications) > 1:
                return find_lca_taxa(all_split_taxa).strip("; ")

            # get return joined string if end is most common
            elif most_common_rank_classifications[0] == "":
                return "; ".join(all_split_taxa[0]).strip("; ")

            else:    
                level_idx += 1
        
        else:
            return ""



# assign_contig_taxonomy()
#
def assign_contig_taxonomy(contig_prots_map, gene_tax_map, gene_bitscore_map, method):
    contig_tax_map = dict()

    for contig, prots in contig_prots_map.items():
        if method == "top":
            best_bitscore = 0
            best_prot = ""
            for prot in prots:
                bitscore = gene_bitscore_map[prot]
                if bitscore > best_bitscore:
                    best_bitscore = bitscore
                    best_prot = prot
            
            best_taxonomy = gene_tax_map[best_prot]

        else:
            taxonomies = [gene_tax_map[prot] for prot in prots]

            all_split_taxa = list()
            for taxa in taxonomies:
                if taxa:
                    taxa += "; "
                    all_split_taxa.append(taxa.split("; "))
                else:
                    all_split_taxa.append("")
            
            if method == "majority":
                best_taxonomy = decide_most_common_taxonomy(all_split_taxa)

            elif method == "lca":
                best_taxonomy = find_lca_taxa(all_split_taxa).strip("; ")

        contig_tax_map[contig] = best_taxonomy
    
    return contig_tax_map
        

# main()
#
def main() -> int:
    args = getargs()

    domtbl_df = read_domtble_to_df(args.hmmer_method, args.domain_table)

    fam_tax_map = read_tax_map(args.fam_tax_map)

    tc_map = read_tc_map(args.tc_map)

    fam_host_map = read_fam_host_file(args.fam_host_file)

    gene_fam_map, gene_bitscore_map = assign_gene_taxonomy(domtbl_df, tc_map)

    gene_tax_map = {gene: fam_tax_map[fam] for gene, fam in gene_fam_map.items()}

    gene_tax_outfile = os.path.join(args.outdir, "gene_tax_assignment.tsv")
    write_gene_tax_to_file(gene_tax_outfile, gene_tax_map, gene_bitscore_map)

    contig_prots_map = get_contigs(list(gene_tax_map.keys()))

    contig_host_counts, contig_host_map = get_contig_host_counts(contig_prots_map, gene_fam_map, fam_host_map)

    contig_host_counts_outfile = os.path.join(args.outdir, "contig_host_counts.tsv")
    write_contig_host_counts(contig_host_counts_outfile, contig_host_counts, contig_host_map)

    for method in ["top", "majority", "lca"]:
        contig_tax_map = assign_contig_taxonomy(contig_prots_map, gene_tax_map, gene_bitscore_map, method)

        contig_tax_outfile = os.path.join(args.outdir, f"{method}_contig_predictions.tsv")
        write_contig_tax_to_file(contig_tax_outfile, contig_tax_map, contig_host_map)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())