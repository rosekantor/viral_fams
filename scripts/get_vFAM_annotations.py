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
import pandas as pd
from collections import defaultdict


# software run
#
#HOME_DIR = '~chivian1'
#HOME_DIR = '/g/g16/chivian1'
#SOFTWARE_DIR = os.path.join (HOME_DIR, 'work', 'software')
#MUSCLE_bin = os.path.join (SOFTWARE_DIR, 'MUSCLE_5.3', 'muscle-linux-x86.v5.3')
#GBLOCKS_bin = os.path.join (SOFTWARE_DIR, 'Gblocks_0.91b', 'Gblocks')
#FASTTREE_bin = os.path.join (SOFTWARE_DIR, 'FastTree_2.1.11', 'FastTreeDbl-linux-x86.v2.1.11')


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--in-sequence-dir", help="")
    parser.add_argument("-s", "--suffix", type=str, help="")
    parser.add_argument("-o", "--outfile", help="")
    
    args = parser.parse_args()
    args_pass = True

    if args.in_sequence_dir is None:
        print ("must specify --{}\n".format('in-sequence-dir'))
        args_pass = False
    if not os.path.exists(args.in_sequence_dir) or not os.path.isdir(args.in_sequence_dir):
        print ("--{} {} must exist and not be empty\n".format('in-sequence-dir', args.in_sequence_dir))
        args_pass = False
    if args.suffix is None:
        print ("must specify --{}\n".format('suffix'))
        args_pass = False
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args

# count_annotations(fasta_file)
#
def count_annotations(fasta_file):
    annotation_counts = defaultdict(int)
    with open(fasta_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith(">"):
                annotation = ""
                getting_annot = False
                for char in line:
                    if char == "[":
                        getting_annot = True
                    elif char == "]":
                        getting_annot = False
                    elif getting_annot:
                        annotation += char
                annotation_counts[annotation] += 1
    
    return annotation_counts

# get_vFam_annotations(in_dir, fasta_suffix)
#
def get_vFam_annotations(in_dir, fasta_suffix):
    fam_annot_counts = defaultdict(dict)
    for db_root, db_dirs, db_files in os.walk(in_dir):
        for fam in db_dirs:
            for fam_root, fam_dirs, fam_files in os.walk(os.path.join(db_root, fam)):
                for file in fam_files:
                    if file.endswith(fasta_suffix):
                        fam_name = file.replace(f"{fasta_suffix}", "")
                        annotation_counts = count_annotations(os.path.join(fam_root, file))
                        fam_annot_counts[fam_name] = dict(sorted(annotation_counts.items(), key=lambda item: item[1], reverse=True))
    return fam_annot_counts

# write_annotations_map_to_file(fam_annotations_map)
#
def write_annotations_map_to_file(outfile, fam_annotations_map):
    with open(outfile, "w") as file:
        for fam, annot_counts in fam_annotations_map.items():
            file.write(f"{fam}\t")
            total = sum(list(annot_counts.values()))
            for annot, count in annot_counts.items():
                file.write(f"{annot} ({count}), ")
            file.write("\n")
            

# main()
#
def main() -> int:
    args = getargs()

    fam_annotations_map = get_vFam_annotations(args.in_sequence_dir, args.suffix)

    write_annotations_map_to_file(args.outfile, fam_annotations_map)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())