#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import random
from tslearn.metrics import dtw

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
    parser = argparse.ArgumentParser(description="Compare novel genome architecture to ICTV genome architectures with similar classifications.")

    parser.add_argument("-p", "--prediction-table", help="final summary report output from parse_hmmscan_table.py")
    parser.add_argument("-g", "--gff-file", help="gff file for input sequences")
    parser.add_argument("-c", "--contig-fams", help="contig fams output from parse_hmmscan_table.py")
    parser.add_argument("-m", "--magnitude-table", help="magnitudes output from get_ref_seqs_architecture-ICTV.py")
    parser.add_argument("-d", "--distance-table", help="raw structure similarity output from get_ref_seqs_architecture-ICTV.py")
    parser.add_argument("-t", "--taxa-fam-sets", help="taxa fams from get_ref_seqs_architecture-ICTV.py")
    parser.add_argument("-o", "--outdir", help="Output directory")
    
    args = parser.parse_args()
    args_pass = True

    if args.prediction_table is None:
        print ("must specify --{}\n".format('prediction-table'))
        args_pass = False
    if not os.path.isfile(args.prediction_table):
        print ("--{} {} must exist and not be empty\n".format('prediction-table', args.prediction_table))
        args_pass = False

    if args.gff_file is None:
        print ("must specify --{}\n".format('gff-file'))
        args_pass = False
    if not os.path.isfile(args.gff_file):
        print ("--{} {} must exist and not be empty\n".format('gff-file', args.gff_file))
        args_pass = False
    
    if args.contig_fams is None:
        print ("must specify --{}\n".format('contig-fams'))
        args_pass = False
    if not os.path.isfile(args.contig_fams):
        print ("--{} {} must exist and not be empty\n".format('contig-fams', args.contig_fams))
        args_pass = False

    if args.magnitude_table is None:
        print ("must specify --{}\n".format('magnitude-table'))
        args_pass = False
    if not os.path.isfile(args.magnitude_table):
        print ("--{} {} must exist and not be empty\n".format('magnitude-table', args.magnitude_table))
        args_pass = False

    if args.distance_table is None:
        print ("must specify --{}\n".format('distance-table'))
        args_pass = False
    if not os.path.isfile(args.distance_table):
        print ("--{} {} must exist and not be empty\n".format('distance-table', args.distance_table))
        args_pass = False

    if args.taxa_fam_sets is None:
        print ("must specify --{}\n".format('taxa-fam-sets'))
        args_pass = False
    if not os.path.isfile(args.taxa_fam_sets):
        print ("--{} {} must exist and not be empty\n".format('taxa-fam-sets', args.taxa_fam_sets))
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


def read_prediction_table(filepath, tax_column):
    predictions = dict()

    df = pd.read_csv(filepath, usecols=["contig", tax_column], sep="\t")

    for i, row in df.iterrows():
        contig = row["contig"]
        taxa = row[tax_column]
        predictions[contig] = taxa
    
    return predictions

# extract_id(attributes_string)
#
def extract_id(attributes_string):
    attributes = attributes_string.split(";")
    for attribute in attributes:
        if attribute.startswith("ID="):
            return attribute.split("=")[1].split("_")[-1]

def read_gff_magnitudes(filepath, included_contigs):
    # contig_id_positions = defaultdict(lambda: defaultdict(list))
    contig_query_magnitudes = defaultdict(lambda: defaultdict(list))

    with open(filepath) as gff_file:
        lines = gff_file.readlines()
        for line in lines:
            if not line.startswith("#"):
                line = line.strip()
                features = line.split("\t")

                if features[2] == "CDS":
                    contig = features[0]
                    if contig in included_contigs:
                        start = int(features[3]) - 1
                        stop = int(features[4])

                        id = extract_id(features[8])
                        query = f"{contig}_{id}"

                        strand = features[6]
                        
                        if strand == "+":
                            magnitude = stop - start
                        elif strand == "-":
                            magnitude = start - stop

                        # contig_id_positions[contig][id].append((start, stop))
                        contig_query_magnitudes[contig][query].append(magnitude)           
    
    '''
    # contig_positions = defaultdict(list)
    contig_magnitudes = defaultdict(list)

    
    for contig, id_positions in contig_id_positions.items():
        for id, positions in id_positions.items():
            contig_positions[contig].append(positions)

    for contig, query_magnitudes in contig_query_magnitudes.items():
        for query, magnitudes in query_magnitudes.items():
            contig_magnitudes[contig].append(magnitudes)
    '''

    return contig_query_magnitudes


def read_contig_fams(filepath):
    contig_query_fams = defaultdict(dict)
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            contig = info[0]
            query = info[1]
            fam = info[2]

            contig_query_fams[contig][query] = fam
    
    return contig_query_fams


def read_magnitude_table(filepath):
    predictions = dict()
    genome_magnitudes = defaultdict(list)
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            tax = info[0]
            genome = info[1]
            magnitudes_info = info[2:]
        
            predictions[genome] = tax
            
            for magnitude in magnitudes_info:
                magnitude = magnitude.strip("[]").split(", ")
                genome_magnitudes[genome].append(magnitude)

    return predictions, genome_magnitudes


def read_distance_table(filepath):
    tax_distances = dict()
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            tax = info[0]
            distance = float(info[1])

            tax_distances[tax] = distance

    return tax_distances


def read_taxa_fam_sets(filepath):
    tax_fams_map = defaultdict(set)
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            tax = info[0]
            fams = {fam for fam in info[1:]}

            tax_fams_map[tax] = fams
    
    return tax_fams_map

def get_taxa_level_groups(tax_map):
    temp_tax_map = {genome: tax.split("; ") for genome, tax in tax_map.items()}

    all_taxa_level_groups = defaultdict(list)

    for genome, taxa_list in temp_tax_map.items():
        for level in range(1, len(taxa_list) + 1):
            all_taxa_level_groups["; ".join(taxa_list[:level])].append(genome)

    return dict(sorted(all_taxa_level_groups.items(), reverse=True))


def get_novel_tax_distances(ictv_predictions, genome_magnitudes, novel_predictions, contig_query_magnitudes):
    novel_tax_distances = defaultdict(dict)

    contig_magnitudes = defaultdict(list)
    for contig, query_magnitudes in contig_query_magnitudes.items():
        for query, magnitudes in query_magnitudes.items():
            contig_magnitudes[contig].append(magnitudes)

    # pad
    max_len = max([len(sub_magnitude) for magnitudes in list(genome_magnitudes.values())+list(contig_query_magnitudes.values()) for sub_magnitude in magnitudes])
    for magnitudes in genome_magnitudes.values():
        for sub_magnitude in magnitudes:
            while len(sub_magnitude) < max_len:
                sub_magnitude.append(0)
    for magnitudes in contig_magnitudes.values():
        for sub_magnitude in magnitudes:
            while len(sub_magnitude) < max_len:
                sub_magnitude.append(0)
    
    ictv_taxa_level_groups = get_taxa_level_groups(ictv_predictions)
    
    
    for novel_contig, predicted_tax in novel_predictions.items():     
        while predicted_tax not in ictv_taxa_level_groups.keys():
            predicted_tax = "; ".join(predicted_tax.split("; ")[:-1])

        if taxanomic_rank_above(predicted_tax, "Family"):
            reference = np.array(contig_magnitudes[novel_contig])
            distances = [dtw(reference, np.array(genome_magnitudes[genome])) for genome in ictv_taxa_level_groups[predicted_tax]]
            distance = sum(distances) / len(distances)
        
            novel_tax_distances[novel_contig][predicted_tax] = distance
    
        else:
            print(f"Skipped {novel_contig}")

    return novel_tax_distances, ictv_taxa_level_groups


def write_novel_tax_distances(outfile, novel_tax_distances):
    with open(outfile, "w") as file:
        for contig, tax_map in novel_tax_distances.items():
            for tax, distance in tax_map.items():
                file.write(f"{contig}\t{tax}\t{round(distance, 2)}\n")


def taxanomic_rank_above(taxa, rank_below):
    novel_rank = identify_rank(taxa)
    return ORDERED_TAXA_RANKS.index(rank_below) <= ORDERED_TAXA_RANKS.index(novel_rank)


def identify_rank(full_taxonomic_string):
    taxonomy_name = full_taxonomic_string.split("; ")[-1]
    if " " in taxonomy_name:
        return "Species"
    for rank, suffixes in RANK_SUFFIX_MAP.items():
        for suffix in suffixes:
            if taxonomy_name.endswith(suffix):
                return rank
    return None


def compare_tax_distances(novel_tax_distances, original_tax_distances):
    good_predictions = set()
    bad_predictions = set()

    for contig, tax_map in novel_tax_distances.items():
        for tax, distance in tax_map.items():
            average_distance = original_tax_distances[tax]
            if average_distance > 0:
                if distance < (average_distance * 1.25):
                    good_predictions.add(contig)
                else:
                    bad_predictions.add(contig)
            else:
                if distance < 5000:
                    good_predictions.add(contig)
                else:
                    bad_predictions.add(contig)

    return good_predictions, bad_predictions


def compare_tax_fams(contig_query_fams, novel_predictions, tax_fams_map):
    flagged_queries = set()

    for contig, query_fam_map in contig_query_fams.items():
        predicted_tax = novel_predictions[contig]
        for query, fam in query_fam_map.items():
            if fam not in tax_fams_map[predicted_tax]:
                flagged_queries.add(query)
    
    return flagged_queries


def find_bad_structure_and_fam(flagged_fam_queries, bad_structural_contigs, contig_query_magnitudes, genome_magnitudes, ictv_taxa_level_groups, novel_predictions, novel_tax_distances):
    bad_queries = defaultdict(list)

    for contig, query_magnitudes in contig_query_magnitudes.items():
        for check_query in query_magnitudes.keys():
            if contig in bad_structural_contigs and len(query_magnitudes) > 1:
                predicted_tax = novel_predictions[contig]
                while predicted_tax not in ictv_taxa_level_groups.keys():
                    predicted_tax = "; ".join(predicted_tax.split("; ")[:-1])

                contig_magnitudes_removed_query = np.array([magnitudes for query, magnitudes in query_magnitudes.items() if query != check_query])
                distances = [dtw(contig_magnitudes_removed_query, np.array(genome_magnitudes[genome])) for genome in ictv_taxa_level_groups[predicted_tax]]
                distance = sum(distances) / len(distances)

                if distance < novel_tax_distances[contig][predicted_tax]:
                    bad_queries[check_query].append("Architecture")
                    print(check_query)
                    print(predicted_tax)
                    print(distance)

            if check_query in flagged_fam_queries:
                bad_queries[check_query].append("FAM")

    return bad_queries


def write_set_to_file(outfile, this_set):
    with open(outfile, "w") as file:
        for name in this_set:
            file.write(f"{name}\n")

def write_dict_to_file(outfile, this_dict):
    with open(outfile, "w") as file:
        for key, val in this_dict.items():
            file.write(f"{key}\t{val}\n")

# main()
#
def main() -> int:
    args = getargs()

    novel_predictions = read_prediction_table(args.prediction_table, "best-hit_taxonomy")

    contig_query_magnitudes = read_gff_magnitudes(args.gff_file, set(novel_predictions.keys()))

    contig_query_fams = read_contig_fams(args.contig_fams)

    ictv_predictions, genome_magnitudes = read_magnitude_table(args.magnitude_table)

    original_tax_distances = read_distance_table(args.distance_table)

    tax_fams_map = read_taxa_fam_sets(args.taxa_fam_sets)

    # fam structure
    # -------------------------------------------------------------------------
    flagged_fam_queries = compare_tax_fams(contig_query_fams, novel_predictions, tax_fams_map)
    # -------------------------------------------------------------------------

    # magnitude structure
    # -------------------------------------------------------------------------
    novel_tax_distances, ictv_taxa_level_groups = get_novel_tax_distances(ictv_predictions, genome_magnitudes, novel_predictions, contig_query_magnitudes)

    novel_tax_dstances_outfile = os.path.join(args.outdir, "contig_tax_distances.tsv")
    write_novel_tax_distances(novel_tax_dstances_outfile, novel_tax_distances)

    good_structural_contigs, bad_structural_contigs = compare_tax_distances(novel_tax_distances, original_tax_distances)

    good_structural_contigs_outfile = os.path.join(args.outdir, "good_structural_contigs.txt")
    write_set_to_file(good_structural_contigs_outfile, good_structural_contigs)

    bad_structural_contigs_outfile = os.path.join(args.outdir, "bad_structural_contigs.txt")
    write_set_to_file(bad_structural_contigs_outfile, bad_structural_contigs)
    # -------------------------------------------------------------------------

    bad_queries = find_bad_structure_and_fam(flagged_fam_queries, bad_structural_contigs, contig_query_magnitudes, genome_magnitudes, ictv_taxa_level_groups, novel_predictions, novel_tax_distances)

    bad_queries_outfile = os.path.join(args.outdir, "bad_queries.txt")
    write_dict_to_file(bad_queries_outfile, bad_queries)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())