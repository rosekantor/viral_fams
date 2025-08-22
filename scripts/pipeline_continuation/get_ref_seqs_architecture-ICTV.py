#!/usr/bin/python3

import sys
import os
import argparse
import re
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
from Bio import SeqIO
import matplotlib.pyplot as plt
import math
from adjustText import adjust_text
import numpy as np
from tslearn.metrics import dtw
import random
import concurrent.futures
import time

# Note: I was working on implementing multiprocessing and I think it works, but I never got back to it

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
    parser = argparse.ArgumentParser(description="Get ICTV sequences for each fam and genome architecture/fam hits.")

    parser.add_argument("-l", "--lin-hits-table", help="output from get_vFAM_lineage_and_host_range.py using the same HMM db")
    parser.add_argument("-d", "--db-name", help="name of HMM db used (used to filter domtbl filenames from --hmm-hits-dir)")
    parser.add_argument("-a", "--hmm-hits-dir", help="genome hits dir")
    parser.add_argument("-g", "--genomes-dir", help="genome sequences dir")
    parser.add_argument("-c", "--cpus", type=int, default=1, help="Number of cpus to use for multiprocessing")
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

# make_dir(path)
#
def make_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)

# read_lin_hits_to_dict(filepath)
#
def read_lin_hits_to_dict(filepath):
    prot_fam_map = defaultdict(lambda: defaultdict(list))
    df = pd.read_csv(filepath, sep="\t", index_col=0)
    genome_hits_map = {fam: row["GENOME_HITS"] for fam, row in df.iterrows()}

    for fam, genomes in genome_hits_map.items():
        genomes = genomes.split(";")
        for genome in genomes:
            accession_info = genome.split(":")
            genome_id = accession_info[0]
            prot_id = accession_info[1]
            prot_fam_map[genome_id][prot_id].append(fam)

    return prot_fam_map


# get_hit_coords(hmm_hits_dir, db_name, prot_fam_map)
#
def get_hit_coords(hmm_hits_dir, db_name, prot_fam_map):
    prot_coords = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    prot_bitscores = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for annot_root, annot_dirs, annot_files in os.walk(hmm_hits_dir):
        for genome in annot_dirs:
            if genome in prot_fam_map.keys():
                for genome_root, genome_dirs, genome_files in os.walk(os.path.join(annot_root, genome)):
                    for file in genome_files:
                        if db_name in file and "domtbl" in file:
                            with open(os.path.join(genome_root, file), 'r') as file:
                                for qresult in SearchIO.parse(file, 'hmmscan3-domtab'):
                                    for hit in qresult.hits:
                                        for hsp in hit.hsps:
                                            query = qresult.id
                                            target = hit.id
                                            bitscore = hsp.bitscore
                                            env_start = hsp.env_start
                                            env_end = hsp.env_end
                                            if query in prot_fam_map[genome].keys():
                                                if target in prot_fam_map[genome][query]:
                                                    prot_coords[genome][query][target].append((env_start, env_end))
                                                    prot_bitscores[genome][query][target].append(bitscore)
    
    return prot_coords, prot_bitscores


# extract_name(attributes_string)
#
def extract_name(attributes_string):
    attributes = attributes_string.split(";")
    for attribute in attributes:
        if attribute.startswith("Name="):
            return attribute.split("=")[1]

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
        
# get_genome_info(genomes_dir, prot_fam_coords)
#
def get_genome_info(genomes_dir, prot_fam_coords):
    fam_seq_map = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    seq_map = defaultdict(dict)
    tax_map = dict()
    gene_pos_map = defaultdict(lambda: defaultdict(list))
    full_genome_region = dict()
    
    for genome_root, genome_dirs, genome_files in os.walk(genomes_dir):
        for genome in genome_dirs:
            for info_root, info_dirs, info_files in os.walk(os.path.join(genome_root, genome)):
                if any([file.endswith("-protein.faa") for file in info_files]):
                    for file in info_files:
                        filepath = os.path.join(info_root, file)
                        
                        if file.endswith("-protein.faa"):
                            fasta_dict = read_fasta_to_dict(filepath)
                            for prot, seq in fasta_dict.items():
                                seq_map[genome][prot] = seq
                                if prot in prot_fam_coords[genome].keys():
                                    for fam, all_coords in prot_fam_coords[genome][prot].items():
                                        for coords in all_coords:
                                            fam_seq_map[fam][genome][prot] = seq[coords[0]:coords[1]]

                        elif file.endswith("-genomic.gbff"):
                            for record in SeqIO.parse(filepath, "genbank"):
                                taxonomy = "; ".join(record.annotations["taxonomy"]) + f"; {record.annotations["organism"]}"
                                tax_map[genome] = taxonomy

                        elif file.endswith("-genomic.gff"):
                            with open(filepath) as gff_file:
                                lines = gff_file.readlines()
                                for line in lines:
                                    if not line.startswith("#"):
                                        line = line.strip()
                                        features = line.split("\t")

                                        if features[2] == "CDS":
                                            start = int(features[3]) - 1
                                            stop = int(features[4])
                                            strand = features[6]
                                            prot = extract_name(features[8])

                                            if strand == "+":
                                                gene_pos_map[genome][prot].append((start, stop))
                                            elif strand == "-":
                                                gene_pos_map[genome][prot].append((stop, start))

                                        elif features[2] == "region":
                                            start = int(features[3]) - 1
                                            stop = int(features[4])
                                            full_genome_region[genome] = (start, stop)

    return fam_seq_map, seq_map, tax_map, gene_pos_map, full_genome_region


def get_genome_darkness(seq_map, prot_fam_map):
    genome_darkess = defaultdict()

    for genome, prot_map in seq_map.items():
        no_hit_count = 0
        for prot, seq in prot_map.items():
            if prot not in prot_fam_map[genome].keys():
                no_hit_count += 1

        genome_darkess[genome] = no_hit_count / len(prot_map)
    
    return genome_darkess


def write_dict_to_file(outfile, map):
    with open(outfile, "w") as file:
        for key, val in map.items():
            file.write(f"{key}\t{val}\n")


def plot_genome_darkness(genome_darkness_map, outfile):    
    plt.hist(list(genome_darkness_map.values()), bins=10, density=False)
    plt.xlabel("% Darkness")
    plt.ylabel("Frequency")
    plt.title(f"Genome Darkness Frequency")
    plt.savefig(outfile)


# write_fam_sequences_to_faa(outdir, seq_map, tax_map)
#
def write_fam_sequences_to_faa(outdir, seq_map, tax_map):
    for fam, genome_map in seq_map.items():
        with open(os.path.join(outdir, f"{fam}-sliced.faa"), "w") as file:
            for genome, prot_map in genome_map.items():
                for prot, seqs in prot_map.items():
                    for seq in seqs:
                        file.write(f">{genome}:{prot}:{tax_map[genome]}\n")
                        file.write(f"{seq}\n")


def filter_lower_fams(prot_bitscores, prot_fam_coords):
    filtered_prot_fam_coords = prot_fam_coords.copy()

    for genome, prot_map in prot_bitscores.items():
        for prot, fam_map in prot_map.items():
            best_fam = None
            highest_bitscore = -math.inf
            for fam, scores in fam_map.items():
                avg_score = sum(scores) / len(scores)
                if avg_score > highest_bitscore:
                    best_fam = fam
                    highest_bitscore = avg_score
            
            for fam in list(filtered_prot_fam_coords[genome][prot].keys()):
                if fam != best_fam:
                    filtered_prot_fam_coords[genome][prot].pop(fam)
    
    return filtered_prot_fam_coords


def write_hit_taxonomy(outfile, prot_fam_coords, tax_map):
    with open(outfile, "w") as file:
        file.write("QUERY\tFAM\tLINEAGE\n")
        for genome, prot_fam_map in prot_fam_coords.items():
            for prot, fam_coords_map in prot_fam_map.items():
                for fam in fam_coords_map.keys():
                    file.write(f"{genome}_{prot}\t{fam}\t{tax_map[genome]}\n")


# get_taxa_level_groups
#
def get_taxa_level_groups(tax_map):
    temp_tax_map = {genome: tax.split("; ") for genome, tax in tax_map.items()}

    all_taxa_level_groups = defaultdict(list)

    for genome, taxa_list in temp_tax_map.items():
        for level in range(1, len(taxa_list)):
            all_taxa_level_groups["; ".join(taxa_list[:level])].append(genome)

    return dict(sorted(all_taxa_level_groups.items(), reverse=True))


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def move_fam_position(fam_position, gap):
    if fam_position[0] <= gap[0]:
        return [(fam_position[0], gap[0]), (gap[1], gap[1] + (fam_position[1] - gap[0]))]
    else:
        return [(fam_position[0] + (gap[1] - gap[0]), fam_position[1] + (gap[1] - gap[0]))]


# get_fam_architecture(prot_fam_coords, gene_pos_map
#
def get_fam_architecture(prot_fam_coords, gene_pos_map):
    new_prot_fam_coords = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for genome, pos_map in gene_pos_map.items():
        for gene, gene_positions in pos_map.items():
            # get the smallest value
            gene_left_end = math.inf
            for coords in gene_positions:
                if coords[0] < gene_left_end:
                    gene_left_end = coords[0]
                if coords[1] < gene_left_end:
                    gene_left_end = coords[1]

            gaps = list()
            order_gene_positions = sorted([(min(start, stop), max(start, stop)) for start, stop in gene_positions])
            if len(order_gene_positions) > 1:
                gap = (order_gene_positions[0][1], )
                for coords in order_gene_positions[1:]:
                    gap += (coords[0], )
                    gaps.append(gap)
                    gap = (coords[1], )

            for fam, all_fam_coords in prot_fam_coords[genome][gene].items():

                for fam_coords in all_fam_coords:
                    # if prot_bitscores[genome][gene][fam] > fam_lowest_thresh_map[fam]:

                    fam_start = gene_left_end + (fam_coords[0] * 3)
                    fam_stop = gene_left_end + (fam_coords[1] * 3)

                    fam_positions = [(fam_start, fam_stop)]

                    # check if it overlaps a gap
                    for gap in gaps:
                        for fam_position in fam_positions:
                            if get_overlap(fam_position, gap) > 0:    
                                fam_positions.remove(fam_position)
                                fam_positions += move_fam_position(fam_position, gap)
                                break
                    
                    for fam_position in fam_positions:
                        new_prot_fam_coords[genome][gene][fam].append(fam_position)

    return new_prot_fam_coords


# write_genome_fam_architecture(outdir, all_taxa_level_groups, new_prot_fam_coords)
#
def write_genome_fam_architecture(outdir, all_taxa_level_groups, new_prot_fam_coords):
    seen_genomes = set()
    with open(os.path.join(outdir, f"fam_architecture.tsv"), "w") as file:
        for tax, genomes in all_taxa_level_groups.items():
            for genome in genomes:
                if genome not in seen_genomes:
                    for gene, fam_map in new_prot_fam_coords[genome].items():
                        for fam, positions in fam_map.items():
                            file.write(f"{tax}\t{genome}\t{gene}\t{fam}")
                            for position in positions:
                                file.write(f"\t{position}")
                            file.write(f"\n")

                    seen_genomes.add(genome)


# write_taxa_fam_sets(outdir, all_taxa_level_groups, prot_fam_map)
#
def write_taxa_fam_sets(outdir, all_taxa_level_groups, prot_fam_map):
    tax_fams = defaultdict(set)
    
    for tax, genomes in all_taxa_level_groups.items():
        for genome in genomes:
            for gene, fams in prot_fam_map[genome].items():
                for fam in fams:
                    tax_fams[tax].add(fam)

    with open(os.path.join(outdir, f"taxa_fams.tsv"), "w") as file:
        for tax, fams in tax_fams.items():
            file.write(f"{tax}")
            order_fams = sorted(list(fams))
            for fam in order_fams:
                file.write(f"\t{fam}")
            file.write(f"\n")


# get_genome_magnitudes(gene_pos_map)
#
def get_genome_magnitudes(gene_pos_map):
    genome_stucture_map = defaultdict(list)

    # get magnitudes
    for genome, pos_map in gene_pos_map.items():
        for i, (gene, gene_positions) in enumerate(pos_map.items()):
            genome_stucture_map[genome].append(list())
            for coords in gene_positions:
                genome_stucture_map[genome][i].append(coords[1] - coords[0])
    
    return genome_stucture_map


# write_genome_magnitudes(outdir, all_taxa_level_groups, genome_stucture_map)
#
def write_genome_magnitudes(outdir, all_taxa_level_groups, genome_stucture_map):
    seen_genomes = set()
    with open(os.path.join(outdir, f"maginitudes.tsv"), "w") as file:
        for tax, genomes in all_taxa_level_groups.items():
            for genome in genomes:
                if genome not in seen_genomes:
                    file.write(f"{tax}\t{genome}")
                    for magnitude in genome_stucture_map[genome]:
                        file.write(f"\t{magnitude}")
                    file.write(f"\n")

                    seen_genomes.add(genome)


# get_taxanomic_level_structural_similarity(genome_stucture_map, all_taxa_level_groups)
#
def get_taxanomic_level_structural_similarity(genome_stucture_map, all_taxa_level_groups, cpus):
    norm_tax_similarity_map = dict()
    raw_tax_similarity_map = dict()

    # pad arrays
    max_len = max([len(sub_array) for genome in genome_stucture_map.keys() for sub_array in genome_stucture_map[genome]])
    max_width = max([len(genome_stucture_map[genome]) for genome in genome_stucture_map.keys()])
    for genome in genome_stucture_map.keys():
        for sub_array in genome_stucture_map[genome]:
            while len(sub_array) < max_len:
                sub_array.append(0)

    for tax, genomes in all_taxa_level_groups.items():
        if taxanomic_rank_above(tax, "Family"):
            distance = 0
            if len(genomes) > 1:
                avg_distances = list()

                for ref_genome in genomes:
                    reference_structure = genome_stucture_map[ref_genome]
                    other_structures = [genome_stucture_map[other_genome] for other_genome in genomes]
                    avg_distances.append(
                        parallel_compute_average_distance(reference_structure, other_structures, cpus)
                    )
            
                distance = sum(avg_distances) / len(avg_distances)

            raw_tax_similarity_map[tax] = distance
            # print(f"{tax}: {distance}")
    
    # normalize
    max_dist = max(raw_tax_similarity_map.values())
    min_dist = min(raw_tax_similarity_map.values())
    for tax, distance in raw_tax_similarity_map.items():
        norm_tax_similarity_map[tax] = (distance - min_dist) / (max_dist - min_dist)

    return norm_tax_similarity_map, raw_tax_similarity_map


def single_dtw(ref_structure, other_structure):
    return dtw(np.array(ref_structure), np.array(other_structure))

def parallel_compute_average_distance(ref_structure, other_structures, cpus):
    distances = list()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = [executor.submit(single_dtw,
                                ref_structure,
                                other_structure
                                ) for other_structure in other_structures]

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            distances.append(result)
    
    return sum(distances) / len(distances)


def jaccard_distance(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return 1 - intersection / union if union > 0 else 0

# get_taxanomic_level_fam_similarity(new_prot_fam_coords, all_taxa_level_groups)
#
def get_taxanomic_level_fam_similarity(new_prot_fam_coords, all_taxa_level_groups):
    tax_similarity_map = dict()

    genome_fam_map = defaultdict(set)
    for genome, prot_fam_map in new_prot_fam_coords.items():
        for gene, fam_map in prot_fam_map.items():
            for fam in fam_map.keys():
                genome_fam_map[genome].add(fam)

    for tax, genomes in all_taxa_level_groups.items():
        distances = list()
        for ref_genome in genomes:
            distances += [jaccard_distance(genome_fam_map[ref_genome], genome_fam_map[other_genome]) for other_genome in genomes]

        distance = sum(distances) / len(distances)
        
        tax_similarity_map[tax] = distance
        # print(f"{tax}: {distance}")

    return tax_similarity_map


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


def plot_level_similarity(tax_structure_similarity_map, tax_fam_similarity_map, outfile):
    level_structure_similarity_map = defaultdict(list)
    level_fam_similarity_map = defaultdict(list)

    for tax in tax_structure_similarity_map.keys():
        struct_distance = tax_structure_similarity_map[tax]
        fam_distance = tax_fam_similarity_map[tax]

        rank = identify_rank(tax)
        level_structure_similarity_map[rank].append(struct_distance)
        level_fam_similarity_map[rank].append(fam_distance)

    avg_level_structure_similarity = list()
    avg_level_fam_similarity = list()
    for label in ORDERED_TAXA_RANKS:
        if label in level_structure_similarity_map.keys():
            avg_level_structure_similarity.append( 1 - (sum(level_structure_similarity_map[label]) / len(level_structure_similarity_map[label])) )
        else:
            avg_level_structure_similarity.append(0)

        if label in level_fam_similarity_map.keys():
            avg_level_fam_similarity.append( 1 - (sum(level_fam_similarity_map[label]) / len(level_fam_similarity_map[label])) )
        else:
            avg_level_fam_similarity.append(0)

    x = np.arange(len(ORDERED_TAXA_RANKS))
    width = 0.35 
    fig, ax = plt.subplots(figsize=(20,15))
    rects1 = ax.bar(x - width/2, avg_level_structure_similarity, width, label='structural similarity')
    rects2 = ax.bar(x + width/2, avg_level_fam_similarity, width, label='fam composition similarity')
    ax.set_ylabel('%')
    ax.set_title('Level Similarity')
    ax.set_xticks(x)
    ax.set_xticklabels(ORDERED_TAXA_RANKS)
    ax.legend()

    plt.savefig(outfile)


# write_architecture_to_file(outdir, prot_fam_coords, gene_pos_map, full_genome_region, tax_map)
#
def write_architecture_to_file(outdir, prot_fam_coords, gene_pos_map, full_genome_region, tax_map):
    
    # for tax, genomes in all_taxa_level_groups.items():
        # if tax == "Viruses; Riboviria; Orthornavirae; Pisuviricota; Stelpaviricetes":


    for genome in full_genome_region.keys():
        tax = tax_map[genome]

        plt.figure(figsize=(15, 10))
        plt.title(f"{tax}")

        texts = list()
        
        pos = 0

        pos += 10

        plt.plot([full_genome_region[genome][0], full_genome_region[genome][1]], [pos, pos], color='black')

        plt.text(full_genome_region[genome][1] + 50, pos, genome, fontsize=12, color='black', ha='center')

        width = full_genome_region[genome][1] - full_genome_region[genome][0]
        height = 0
        for gene, gene_positions in gene_pos_map[genome].items():
            height += 10 * ((len(prot_fam_coords[genome][gene]))+1)

        for gene, gene_positions in gene_pos_map[genome].items():
            pos += 10

            for coords in gene_positions:
                gene_section_start = coords[0]
                gene_section_stop = coords[1]
                plt.arrow(
                    gene_section_start, 
                    pos, 
                    gene_section_stop - gene_section_start, 
                    0, 
                    head_width=height * 0.01, 
                    head_length=width * 0.01,
                    length_includes_head = True, 
                    color='red'
                )
                
            # plt.text(gene_section_stop + 50, pos, gene, fontsize=12, color='black', ha='center')

            for fam, all_fam_coords in prot_fam_coords[genome][gene].items():
                pos += 10

                for fam_coords in all_fam_coords:
                    plt.arrow(
                        fam_coords[0], 
                        pos, 
                        fam_coords[1] - fam_coords[0], 
                        0, 
                        # head_width=height * 0.005, 
                        # head_length=width * 0.005, 
                        length_includes_head = True, 
                        color='blue'
                    )

                texts.append(plt.text(fam_coords[1], pos, fam, fontsize=10, color='black', ha='left'))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))
        plt.savefig(os.path.join(outdir, f"{tax.replace(" ", "_").replace("/", "_")};_{genome}.png"))
        # plt.show()

        plt.close()


# main()
#
def main() -> int:
    t_start = time.time()

    args = getargs()

    # genome -> protein -> fams_list
    prot_fam_map = read_lin_hits_to_dict(args.lin_hits_table)

    # genome -> protein -> fam: [(start1, end1), ...], genome -> protein -> fam: [bitscore1, ...]
    prot_fam_coords, prot_bitscores = get_hit_coords(args.hmm_hits_dir, args.db_name, prot_fam_map)

    (   
        fam_seq_map,        # fam -> genome -> prot -> [seq1, ...]
        seq_map,            # genome -> prot -> seq
        tax_map,            # genome -> tax
        gene_pos_map,     # genome -> prot -> [(start_1, end_1), ...]
        full_genome_region  # genome -> (start, end)
    ) = get_genome_info(args.genomes_dir, prot_fam_coords)

    genome_darkness_map = get_genome_darkness(seq_map, prot_fam_map)

    darkness_outfile = os.path.join(args.outdir, "genome_darkness.tsv")
    write_dict_to_file(darkness_outfile, genome_darkness_map)

    outfile = os.path.join(args.outdir, "genome_darkness.png")
    plot_genome_darkness(genome_darkness_map, outfile)

    fastas_outdir = os.path.join(args.outdir, "sequences")
    make_dir(fastas_outdir)
    write_fam_sequences_to_faa(fastas_outdir, fam_seq_map, tax_map)

    prot_fam_coords = filter_lower_fams(prot_bitscores, prot_fam_coords)

    hits_outfile = os.path.join(args.outdir, "hits.tsv")
    write_hit_taxonomy(hits_outfile, prot_fam_coords, tax_map)

    all_taxa_level_groups = get_taxa_level_groups(tax_map)

    # genome -> protein -> fam: [(start1, end1), ...]
    prot_fam_coords = get_fam_architecture(prot_fam_coords, gene_pos_map)
    # write_genome_fam_architecture(args.outdir, all_taxa_level_groups, prot_fam_coords)
    write_taxa_fam_sets(args.outdir, all_taxa_level_groups, prot_fam_map)

    genome_stucture_map = get_genome_magnitudes(gene_pos_map)
    write_genome_magnitudes(args.outdir, all_taxa_level_groups, genome_stucture_map)

    norm_tax_structure_similarity_map, raw_tax_structure_similarity_map = get_taxanomic_level_structural_similarity(genome_stucture_map, all_taxa_level_groups, args.cpus)

    outfile = os.path.join(args.outdir, "raw_taxonomic_level_structure_similarity.tsv")
    write_dict_to_file(outfile, raw_tax_structure_similarity_map)

    outfile = os.path.join(args.outdir, "normalized_taxonomic_level_structure_similarity.tsv")
    write_dict_to_file(outfile, norm_tax_structure_similarity_map)

    tax_fam_similarity_map = get_taxanomic_level_fam_similarity(prot_fam_coords, all_taxa_level_groups)

    outfile = os.path.join(args.outdir, "taxonomic_level_fam_similarity.tsv")
    write_dict_to_file(outfile, tax_fam_similarity_map)

    outfile = os.path.join(args.outdir, "level_similarity.png")
    plot_level_similarity(norm_tax_structure_similarity_map, tax_fam_similarity_map, outfile)

    # write genome architecture to files
    architecture_outdir = os.path.join(args.outdir, "architecture")
    make_dir(architecture_outdir)
    # write_architecture_to_file(architecture_outdir, prot_fam_coords, gene_pos_map, full_genome_region, tax_map) # <-- todo: this gets a filename too long error

    t_end = time.time()
    print ("DONE")
    print(f"Time: {t_end - t_start}")
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

