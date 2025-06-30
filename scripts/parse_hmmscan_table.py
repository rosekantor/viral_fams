#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
import matplotlib.pyplot as plt
import numpy as np


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

METHODS = ["lca", "majority", "best-hit"]

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--domain-table", help="")
    parser.add_argument("-l1", "--lin-thresh-table", help="")
    parser.add_argument("-l2", "--lin-host-table", help="")
    # parser.add_argument("-m", "--contig-taxa-method", default="lca", help="Options: lca, majority, best-hit")
    parser.add_argument("-f", "--fasta", help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.domain_table is None:
        print ("must specify --{}\n".format('domain-table'))
        args_pass = False
    if not os.path.isfile(args.domain_table):
        print ("--{} {} must exist and not be empty\n".format('domain-table', args.domain_table))
        args_pass = False
    if args.lin_thresh_table is None:
        print ("must specify --{}\n".format('lin-thresh-table'))
        args_pass = False
    if not os.path.isfile(args.lin_thresh_table):
        print ("--{} {} must exist and not be empty\n".format('lin-thresh-table', args.lin_thresh_table))
        args_pass = False
    if args.lin_host_table is None:
        print ("must specify --{}\n".format('lin-host-table'))
        args_pass = False
    if not os.path.isfile(args.lin_host_table):
        print ("--{} {} must exist and not be empty\n".format('lin-host-table', args.lin_host_table))
        args_pass = False
    '''
    if args.contig_taxa_method not in ["lca", "majority", "best-hit"]:
        print ("--{} must be a valid method option\n".format('contig-taxa-method'))
        args_pass = False
    '''

    if args.fasta is None:
        print ("must specify --{}\n".format('fasta'))
        args_pass = False
    if not os.path.isfile(args.fasta):
        print ("--{} {} must exist and not be empty\n".format('fasta', args.fasta))
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
def read_domtble_to_df(filepath):
    hits = list()

    with open(filepath, 'r') as file:
        for qresult in SearchIO.parse(file, 'hmmscan3-domtab'):
            for hit in qresult.hits:
                for hsp in hit.hsps:
                    query = qresult.id
                    target = hit.id
                    bitscore = hsp.bitscore
                    evalue = hsp.evalue
                    env_start = hsp.env_start
                    env_end = hsp.env_end

                    hits.append((query, target, bitscore, evalue, env_start, env_end))

    return pd.DataFrame(hits, columns=["query name", "target name", "score", "evalue", "env start", "env end"])

# read_lin_thresh_to_dict(filepath)
#
def read_lin_thresh_to_dict(filepath):
    df = pd.read_csv(filepath, sep="\t", usecols=[0, 1], index_col=0)
    
    # use first thresh for each fam
    first_df = df[~df.index.duplicated(keep='first')]
    last_df = df[~df.index.duplicated(keep='last')]

    return {fam: info["BITSCORE_THRESH"] for fam, info in first_df.iterrows()}, {fam: info["BITSCORE_THRESH"] for fam, info in last_df.iterrows()}

# filter_hit_table(domtbl_df, fam_thresh_map)
#
def filter_hit_table(domtbl_df, fam_thresh_map, invalid_outfile):
    query_fam_map = dict()
    query_bitscore_map = dict()
    query_coords_map = defaultdict(lambda: {"env start": 0, "env end": 0})
    target_not_found_list = list()
    evalue_thresh = 0.001

    query_grouped_tbl_df = domtbl_df.groupby("query name")
    for query, df in query_grouped_tbl_df:
        target_found = False
        df.sort_values(by="score", ascending=False)
        for i, info in df.iterrows():
            target = info["target name"]
            bitscore = info["score"]
            evalue = info["evalue"]
            env_start = info["env start"]   # Note: this automatically subtracts one from raw data
            env_end = info["env end"]
            if evalue < evalue_thresh and target in fam_thresh_map.keys():
                if bitscore > fam_thresh_map[target]:
                    query_fam_map[query] = target
                    query_bitscore_map[query] = bitscore
                    query_coords_map[query]["env start"] = env_start
                    query_coords_map[query]["env end"] = env_end

                    target_found = True
                    break
        
        if not target_found:
            target_not_found_list.append(query)
    
    with open(invalid_outfile, "w") as outfile:
        for invalid_query in target_not_found_list:
            outfile.write(f"{invalid_query}\n")
    
    return query_fam_map, query_bitscore_map, query_coords_map

# read_lin_hosts_to_dicts(filepath)
#
def read_lin_hosts_to_dicts(filepath):
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    df["bacteria+phage"] = df["bacteria"] + df["phage"]
    df.drop(["bacteria", "phage"], axis=1, inplace=True)

    taxonomy_map = {fam: row["HITS_LCA_NODE"] for fam, row in df.iterrows()}
    
    hosts = df.columns[2:]
    host_count_map = defaultdict(lambda: {host: 0 for host in hosts})
    for fam, row in df.iterrows():
        for host in hosts:
            host_count_map[fam][host] += row[host]

    return taxonomy_map, host_count_map

# connect_query_to_fam_map(query_hits_map, fam_info_map)
#
def connect_query_to_fam_map(query_hits_map, fam_info_map):
    query_hosts_map = dict()
    for query, fam in query_hits_map.items():
        query_hosts_map[query] = fam_info_map[fam]
    
    return query_hosts_map

# write_query_info_to_tsv(outfile, query_fam_map, query_bitscore_map, query_taxa_map, query_hosts_map)
#
def write_query_info_to_tsv(outfile, query_fam_map, query_bitscore_map, query_taxa_map, query_hosts_map):
    out_data = list()
    column_names = ["QUERY", "FAM", "BITSCORE", "LINEAGE"] + [host for host in query_hosts_map[list(query_hosts_map.keys())[0]].keys()]

    for query, fam in query_fam_map.items():
        host_data = list()
        for count in query_hosts_map[query].values():
            host_data.append(count)
        
        out_data.append((query, fam, query_bitscore_map[query], query_taxa_map[query])+tuple(host_data))
    
    pd.DataFrame(out_data, columns=column_names).to_csv(outfile, sep="\t", index=False)

# get_contig_id(query_name)
#
def get_contig_id(query_name):
    return "_".join(query_name.split("_")[:-1])

# find_lca_taxa(all_taxa)
#
def find_lca_taxa(all_taxa):
    all_split_taxa = [taxa.split("; ") for taxa in all_taxa]
    lca_index = 0
    while True:
        if not all(lca_index < len(taxa) for taxa in all_split_taxa):
            break
        if not all(taxa[lca_index] == all_split_taxa[0][lca_index] for taxa in all_split_taxa):
            break
        lca_index += 1

    return "; ".join(all_split_taxa[0][:lca_index])

# find_majority_taxa(all_taxa)
#
def find_majority_taxa(all_taxa):
    taxa_counts = defaultdict(int)
    for taxa in all_taxa:
        taxa_counts[taxa] += 1

    greatest_count = 0
    majority_taxa = ""
    for taxa, count in taxa_counts.items():
        if count > greatest_count:
            greatest_count = count
            majority_taxa = taxa

    return majority_taxa

# find_best_hit_taxa(queries, query_bitscore_map, query_taxa_map)
#
def find_best_hit_taxa(queries, query_bitscore_map, query_taxa_map):
    best_bitscore = 0
    best_query = ""
    for query in queries:
        query_bitscore = query_bitscore_map[query]
        if query_bitscore > best_bitscore:
            best_bitscore = query_bitscore
            best_query = query

    return query_taxa_map[best_query]

# determine_contig_taxa(query_taxa_map)
#
def determine_contig_taxa(method, query_taxa_map, contig_all_queries_map, query_bitscore_map):
    contig_all_taxa_map = defaultdict(list)
    for query, taxa in query_taxa_map.items():
        contig_id = get_contig_id(query)
        contig_all_taxa_map[contig_id].append(taxa)
        
    contig_taxa_map = dict()
    for contig, all_taxa in contig_all_taxa_map.items():
        if method == "lca":
            contig_taxa_map[contig] = find_lca_taxa(all_taxa)
        elif method == "majority":
            contig_taxa_map[contig] = find_majority_taxa(all_taxa)
        elif method == "best-hit":
            contig_taxa_map[contig] = find_best_hit_taxa(contig_all_queries_map[contig], query_bitscore_map, query_taxa_map)
        else:
            raise NameError("invalid contig taxa grouping method")
        
    return contig_taxa_map

# write_contig_taxa_to_tsv(outfile, contig_taxa_map)
#
def write_contig_taxa_to_tsv(outfile, contig_taxa_map):
    out_data = [(contig, taxa) for contig, taxa in contig_taxa_map.items()]
    pd.DataFrame(out_data, columns=["contig", "taxonomy"]).to_csv(outfile, sep="\t", index=False)

# identify_rank(taxonomy_name)
#
def identify_rank(taxonomy_name):
    if " " in taxonomy_name:
        return "Species"
    for rank, suffixes in RANK_SUFFIX_MAP.items():
        for suffix in suffixes:
            if taxonomy_name.endswith(suffix):
                return rank
    return None

# create_taxa_summary(summary_outfile, classification_outfile, contig_taxa_map)
#
def create_taxa_summary(summary_outfile, classification_outfile, contig_taxa_map):
    taxa_rank_counts = defaultdict(int)
    taxa_sum_map = dict()
    taxa_lists_map = defaultdict(list)

    for name, taxa in contig_taxa_map.items():
        split_taxa = taxa.split("; ")
        for taxa_name in split_taxa:
            rank = identify_rank(taxa_name)
            if rank:
                taxa_rank_counts[rank] += 1
                taxa_lists_map[name].append(rank)
            else:
                raise NameError(f"'{taxa_name}' from '{taxa}' does not have a mappable taxonomic rank")
        
    sorted_taxa_level_counts = {rank: taxa_rank_counts[rank] for rank in ORDERED_TAXA_RANKS if rank in taxa_rank_counts.keys()}
    total = len(contig_taxa_map)
    with open(summary_outfile, 'w') as file:
        for rank, count in sorted_taxa_level_counts.items():
            percentage = round((count / total) * 100, 2)
            taxa_sum_map[rank] = percentage
            file.write(f"{rank}: {percentage} %\n")

    write_taxa_classification(classification_outfile, taxa_lists_map)
    
    return taxa_sum_map, taxa_lists_map

# write_taxa_classification(classification_outfile, taxa_lists_map)
#
def write_taxa_classification(classification_outfile, taxa_lists_map):
    outdata = list()
    columns = ["Name"] + list(RANK_SUFFIX_MAP.keys()) + ["Sum"]

    for name, taxa_ranks in taxa_lists_map.items():
        sum = 0
        row_data = [name]
        for rank in RANK_SUFFIX_MAP.keys():
            found_class = 1 if rank in taxa_ranks else 0
            row_data.append(found_class)

            sum += found_class
        
        row_data.append(sum)
        
        outdata.append(tuple(row_data))
    
    pd.DataFrame(outdata, columns=columns).to_csv(classification_outfile, sep="\t", index=False)

# create_host_summary(outfile, contig_queries_map, query_hosts_map)
#
def create_host_summary(outfile, contig_queries_map, query_hosts_map):
    contig_host_counts_map = defaultdict(lambda: defaultdict(int))
    contig_host_prediction_map = defaultdict(lambda: defaultdict(float))
    for contig, queries in contig_queries_map.items():
        for query in queries:
            for host, count in query_hosts_map[query].items():
                contig_host_counts_map[contig][host] += count
    
    with open(outfile, 'w') as file:
        for contig, host_map in contig_host_counts_map.items():
            file.write(f"{contig}\t")

            sorted_host_map = dict(sorted(host_map.items(), key=lambda item: item[1], reverse=True))
            total = sum(list(sorted_host_map.values()))
            top = True
            for host, count in sorted_host_map.items():
                if top:
                    contig_host_prediction_map[contig] = host
                    top = False
                if count > 0:
                    percentage = round((count / total) * 100, 2)
                    file.write(f"{host} ({percentage} %), ")

            file.write("\n")
    
    return contig_host_prediction_map

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
                name = line.split(" # ")[0][1:]
                seq = ""
            else:
                seq += line
                
        if name and seq:
            seq_map[name] = seq
    
    return seq_map

# write_dict_to_fasta(filepath, seq_map)
#
def write_dict_to_fasta(filepath, seq_map):
    with open(filepath, "w") as file:
        for name, seq in seq_map.items():
            file.write(f">{name}\n")
            file.write(f"{seq}\n")

# create_fam_fastas(outdir, fasta_filepath, query_fam_map, query_coords_map)
#
def create_fam_fastas(outdir, fasta_filepath, query_fam_map, query_coords_map):
    seq_map = read_fasta_to_dict(fasta_filepath)
    seq_map = {name.split(" ")[0]: seq for name, seq in seq_map.items()}  # rename keys because hmmscan changes sequence names
    
    sliced_query_name_map = dict()
    fam_query_map = defaultdict(list)
    for query, fam in query_fam_map.items():
        fam_query_map[fam].append(query)
    
    for fam, queries in fam_query_map.items():
        fam_out_fasta = os.path.join(outdir, f"{fam}-sliced.faa")

        fam_seq_map = dict()
        for query in queries:
            start = query_coords_map[query]["env start"]
            stop = query_coords_map[query]["env end"]
            sliced_query_name = f"{query}-slice:{start}-{stop}"
            fam_seq_map[sliced_query_name] = seq_map[query][start:stop]
            sliced_query_name_map[query] = sliced_query_name

        write_dict_to_fasta(fam_out_fasta, fam_seq_map)
    
    return sliced_query_name_map

# write_query_bitscore_map(bitscore_map_outdir, sliced_query_name_map, query_bitscore_map)
#
def write_query_bitscore_map(bitscore_map_outdir, sliced_query_name_map, query_bitscore_map):
    with open(bitscore_map_outdir, "w") as file:
        for name, sliced_name in sliced_query_name_map.items():
            file.write(f"{sliced_name}\t{query_bitscore_map[name]}\n")

# plot_taxa_specificity(taxa_sum_vis_out, fam_taxa_sum_map, contig_taxa_sum_map)
#
def plot_taxa_specificity(taxa_sum_vis_out, fam_taxa_sum_map, contig_taxa_sum_map):
    fam_taxa_names = set(fam_taxa_sum_map.keys())
    contig_taxa_names = set(contig_taxa_sum_map.keys())
    labels = list(RANK_SUFFIX_MAP.keys())

    fam_taxa_percentages = list()
    contig_taxa_percentages = list()
    for label in labels:
        if label in fam_taxa_names:
            fam_taxa_percentages.append(fam_taxa_sum_map[label])
        else:
            fam_taxa_percentages.append(0)

        if label in contig_taxa_names:
            contig_taxa_percentages.append(contig_taxa_sum_map[label])
        else:
            contig_taxa_percentages.append(0)

    x = np.arange(len(labels))
    width = 0.35 
    fig, ax = plt.subplots(figsize=(20,15))
    rects1 = ax.bar(x - width/2, fam_taxa_percentages, width, label='Fam rank percentage')
    rects2 = ax.bar(x + width/2, contig_taxa_percentages, width, label='Contig rank percentage')
    ax.set_ylabel('%')
    ax.set_title('Taxanomic Classification Percentage')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.savefig(taxa_sum_vis_out)

# get_contigs(queries)
#
def get_contigs(queries):
    contig_all_queries_map = defaultdict(list)
    for query in queries:
        contig_id = get_contig_id(query)
        contig_all_queries_map[contig_id].append(query)
    
    return contig_all_queries_map

# get_contig_confidence_score(contig_queries_map, query_bitscore_map)
#
def get_contig_confidence_score(contig_queries_map, query_bitscore_map):
    contig_confidence_map = dict()

    total_contig_bitscore_map = defaultdict(int)
    for contig, queries in contig_queries_map.items():
        for query in queries:
            total_contig_bitscore_map[contig] += query_bitscore_map[query]
    
    min_bitscore_sum = min(list(total_contig_bitscore_map.values()))
    max_bitscore_sum = max(list(total_contig_bitscore_map.values()))

    for contig, total_score in total_contig_bitscore_map.items():
        confidence = (total_score - min_bitscore_sum) / (max_bitscore_sum - min_bitscore_sum)
        contig_confidence_map[contig] = confidence
    
    return contig_confidence_map

# create_contig_summary_report(contig_summary_report_outfile, method_contig_taxonomy_map, contig_host_map, contig_confidence_map)
#
def create_contig_summary_report(contig_summary_report_outfile, method_contig_taxonomy_map, contig_host_map, contig_confidence_map):
    out_data = list()

    for contig, host_percentages in contig_host_map.items():
        row_data = [contig]
        for method in METHODS:
            row_data.append(method_contig_taxonomy_map[method][contig])

        out_data.append(tuple(row_data + [host_percentages, round(contig_confidence_map[contig], 2)]))
    
    pd.DataFrame(out_data, columns=["contig"]+[f"{method}_taxonomy" for method in METHODS]+["host_prediction", "confidence"]).to_csv(contig_summary_report_outfile, sep="\t", index=False)

# make_dir(path)
#
def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)

# main()
#
def main() -> int:
    args = getargs()

    domtbl_df = read_domtble_to_df(args.domain_table)
    # domtbl_df.to_csv("/p/lustre1/golez1/tmp_hmm_hits.tsv", sep="\t", index=False)
    # domtbl_df = pd.read_csv("/p/lustre1/golez1/hmm_hits.tsv", sep="\t")
    fam_first_thresh_map, fam_last_thresh_map = read_lin_thresh_to_dict(args.lin_thresh_table)
    fam_taxa_map, fam_host_count_map = read_lin_hosts_to_dicts(args.lin_host_table)

    fam_taxa_sum_outfile = os.path.join(args.outdir, "fam_taxa_summary.txt")
    fam_taxa_classification_outfile = os.path.join(args.outdir, "fam_taxa_classification.tsv")
    fam_taxa_sum_map, taxa_fams_map = create_taxa_summary(fam_taxa_sum_outfile, fam_taxa_classification_outfile, fam_taxa_map)

    for which_thresh in ["first", "second"]:
        thresh_outdir = os.path.join(args.outdir, f"{which_thresh}_bitscore_thresh")
        make_dir(thresh_outdir)

        invalid_queries_outfile = os.path.join(thresh_outdir, "invalid_queries.txt")
        if which_thresh == "first":
            query_fam_map, query_bitscore_map, query_coords_map = filter_hit_table(domtbl_df, fam_first_thresh_map, invalid_queries_outfile)
        elif which_thresh == "second":
            query_fam_map, query_bitscore_map, query_coords_map = filter_hit_table(domtbl_df, fam_last_thresh_map, invalid_queries_outfile)
        query_taxa_map = connect_query_to_fam_map(query_fam_map, fam_taxa_map)
        query_hosts_map = connect_query_to_fam_map(query_fam_map, fam_host_count_map)

        query_info_outfile = os.path.join(thresh_outdir, "query_fam_lineage_host_info.tsv")
        write_query_info_to_tsv(query_info_outfile, query_fam_map, query_bitscore_map, query_taxa_map, query_hosts_map)

        contig_queries_map = get_contigs(list(query_taxa_map.keys()))

        method_contig_taxonomy_map = defaultdict(dict)
        for contig_taxa_method in METHODS:
            taxa_method_outdir = os.path.join(thresh_outdir, f"contig_taxa_info_{contig_taxa_method}")
            make_dir(taxa_method_outdir)
            
            contig_taxa_map = determine_contig_taxa(contig_taxa_method, query_taxa_map, contig_queries_map, query_bitscore_map)

            contig_taxa_outfile = os.path.join(taxa_method_outdir, "contig_taxa.tsv")
            write_contig_taxa_to_tsv(contig_taxa_outfile, contig_taxa_map)
            
            contig_taxa_sum_outfile = os.path.join(taxa_method_outdir, "contig_taxa_summary.txt")
            contig_taxa_classification_outfile = os.path.join(taxa_method_outdir, "contig_taxa_classification.tsv")
            contig_taxa_sum_map, taxa_contig_lists_map = create_taxa_summary(contig_taxa_sum_outfile, contig_taxa_classification_outfile, contig_taxa_map)

            taxa_sum_vis_out = os.path.join(taxa_method_outdir, "taxa_sum_visualization.png")
            plot_taxa_specificity(taxa_sum_vis_out, fam_taxa_sum_map, contig_taxa_sum_map)

            method_contig_taxonomy_map[contig_taxa_method] = contig_taxa_map

        contig_host_sum_outfile = os.path.join(thresh_outdir, "contig_host_summary.txt")
        contig_host_map = create_host_summary(contig_host_sum_outfile, contig_queries_map, query_hosts_map)

        fastas_outdir = os.path.join(thresh_outdir, "sequences")
        make_dir(fastas_outdir)
        sliced_query_name_map = create_fam_fastas(fastas_outdir, args.fasta, query_fam_map, query_coords_map)

        bitscore_map_outdir = os.path.join(thresh_outdir, "bitscore_map.tsv")
        write_query_bitscore_map(bitscore_map_outdir, sliced_query_name_map, query_bitscore_map)

        contig_confidence_map = get_contig_confidence_score(contig_queries_map, query_bitscore_map)

        contig_summary_report_outfile = os.path.join(thresh_outdir, "contig_final_summary_report.tsv")
        create_contig_summary_report(contig_summary_report_outfile, method_contig_taxonomy_map, contig_host_map, contig_confidence_map)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

