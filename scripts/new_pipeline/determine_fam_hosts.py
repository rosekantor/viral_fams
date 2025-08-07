#!/usr/bin/python3

import sys
import os
import argparse
import pandas as pd
from collections import defaultdict
import numpy as np
import requests
from bs4 import BeautifulSoup

# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-m", "--members-file", help="")
    parser.add_argument("-f", "--faa-file", help="")
    parser.add_argument("-s", "--species-file", help="")
    parser.add_argument("-o", "--outdir", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.members_file is None:
        print ("must specify --{}\n".format('members-file'))
        args_pass = False
    elif not os.path.exists(args.members_file) or \
         not os.path.isfile(args.members_file) or \
         not os.path.getsize(args.members_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('members-file', args.members_file))
        args_pass = False
    
    if args.faa_file is None:
        print ("must specify --{}\n".format('faa-file'))
        args_pass = False
    elif not os.path.exists(args.faa_file) or \
         not os.path.isfile(args.faa_file) or \
         not os.path.getsize(args.faa_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('faa-file', args.faa_file))
        args_pass = False

    if args.species_file is None:
        print ("must specify --{}\n".format('species-file'))
        args_pass = False
    elif not os.path.exists(args.species_file) or \
         not os.path.isfile(args.species_file) or \
         not os.path.getsize(args.species_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('species-file', args.species_file))
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


def read_fam_members(filepath):
    members_map = defaultdict(list)

    df = pd.read_csv(filepath, sep="\t", usecols=["#GroupName", "ProteinIDs"])

    for i, row in df.iterrows():
        group = row["#GroupName"]
        protein_ids = row["ProteinIDs"].split(",")

        for i in range(len(protein_ids)):
            if ".sequence" in protein_ids[i]:
                protein_ids[i] = ".".join(protein_ids[i].split(".")[:-1])

        members_map[group] = protein_ids

    return members_map


def read_gene_species(faa_filepath):
    species_map = dict()
    duplicates = defaultdict(int)

    with open(faa_filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip("\n")
            if line.startswith(">"):
                gene = line[1:].split()[0]

                ''' Commenting this out because vogdb is not consistent with duplicate naming convention
                # Note: VOGDB adds .sequence# for duplicate protein names in the members file
                if gene in species_map.keys() or gene in duplicates.keys():
                    if gene not in duplicates.keys():
                        duplicates[gene] = 1
                        species_map[f"{gene}.sequence{duplicates[gene]}"] = species_map.pop(gene)

                    duplicates[gene] += 1
                    gene = f"{gene}.sequence{duplicates[gene]}"
                '''
                
                species = ""
                bracket_excess = 0
                for char in line:
                    if char == "]":
                        bracket_excess -= 1

                    if bracket_excess:
                        species += char

                    if char == "[":
                        if not bracket_excess:
                            species = ""
                        bracket_excess += 1

                species_map[gene] = species
    
    return species_map


def read_species_tax_id(filepath):
    tax_id_map = defaultdict(list)

    df = pd.read_csv(filepath, sep="\t", usecols=["#species name", "taxon id"])

    for i, row in df.iterrows():
        species = row["#species name"]
        tax_id = row["taxon id"]

        tax_id_map[species] = tax_id

    return tax_id_map

# NOTE: This could (and should) be parallelized
def get_species_host(tax_id_map):
    host_map = dict()

    for species, tax_id in tax_id_map.items():
        url = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=" + str(tax_id)
        response = requests.get(url)
        response.raise_for_status() 

        soup = BeautifulSoup(response.text, 'html.parser')

        text = soup.get_text(separator="\n", strip=True)
        split_text = text.split("\n")

        host = "None"
        for i, line in enumerate(split_text):
            if line.startswith("Host: "):
                host = line[len("Host: "):]
                break

        host_map[species] = host

        print(f"{species}: {host}")
    
    return host_map


def write_dict_to_file(outfile, map):
    with open(outfile, "w") as file:
        for key, val in map.items():
            file.write(f"{key}\t{val}\n")


def read_dict(filepath):
    map = dict()
    with open(filepath, "r") as file:
        lines = file.readlines()
        for line in lines:
            info = line.strip().split("\t")
            map[info[0]] = info[1]
    
    return map


def get_fam_species_counts(fam_members_map, gene_species_map, species_host_map):
    fam_species_counts = defaultdict(lambda: defaultdict(int))

    for fam, members in fam_members_map.items():
        for member in members:
            species = gene_species_map[member]
            host = species_host_map[species]

            fam_species_counts[fam][host] += 1

    return fam_species_counts


def write_fam_species_host_counts(outfile, fam_species_counts):
    all_species = sorted(list(set(species for species_map in fam_species_counts.values() for species in species_map.keys())))

    with open(outfile, "w") as file:
        file.write("Group\tBest Species")
        for species in all_species:
            file.write(f"\t{species}")
        file.write("\n")

        for fam, species_counts in fam_species_counts.items():
            best_host = max(species_counts, key=species_counts.get)

            file.write(f"{fam}\t{best_host}")
            for species in all_species:
                if species in species_counts.keys():
                    file.write(f"\t{species_counts[species]}")
                else:
                    file.write(f"\t{0}")
            file.write("\n")


def main():
    args = getargs()

    fam_members_map = read_fam_members(args.members_file)

    gene_species_map = read_gene_species(args.faa_file)

    species_tax_id_map = read_species_tax_id(args.species_file)

    # species_host_map = get_species_host(species_tax_id_map)
    # write_dict_to_file("/p/lustre1/golez1/vogdb_231_species_host_map.tsv", species_host_map)
    species_host_map = read_dict("/p/lustre1/golez1/vogdb_231_species_host_map.tsv")
    
    fam_species_counts = get_fam_species_counts(fam_members_map, gene_species_map, species_host_map)

    fam_species_counts_outfile = os.path.join(args.outdir, "fam_species_counts.tsv")
    write_fam_species_host_counts(fam_species_counts_outfile, fam_species_counts)

    return 0


if __name__ == "__main__":
    main()