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
import re


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="combine tsvs for hmmer scans against viral gene families")

    parser.add_argument("-a", "--annotfile", help="hmmer domtbl annotation file")
    parser.add_argument("-o", "--outtsvfile", help="output merged hits tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit (-1)

    if args.annotfile is None:
        print ("must specify --{}\n".format('annotfile'))
        args_pass = False
    elif not os.path.exists(args.annotfile) or \
         not os.path.isfile(args.annotfile) or \
         not os.path.getsize(args.annotfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('annotfile', args.annotfile))
        args_pass = False
        
    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_genome_ids()
#
def read_genome_ids (genome_ids_file):
    genome_ids = []

    with open (genome_ids_file, 'r') as gid_h:
        for gid_line in gid_h:
            gid_line = gid_line.strip()
            if not gid_line or gid_line.startswith('#'):
                continue
            gid = gid_line.split()[0]
            genome_ids.append(gid)

    return genome_ids


# read_query_hits ()
#
def read_query_hits (annotfile, dom_i_evalue_thresh, fields):
    query_hits = dict()

    hmmer_dom_annot_path = annotfile
    hmmer_dom_annot_buf = get_hmmer_dom_file_buf (hmmer_dom_annot_path)

    hits = []
    in_table = False
    for hd_line in hmmer_dom_annot_buf:

        hd_line = hd_line.strip()
        if hd_line.startswith("# target name"):
            in_table = True
            continue
        elif not hd_line or hd_line.startswith('#'):
            continue
        elif not in_table:
            continue
            
        rec = hd_line.split()
        hit = {}
        for f_i,f in enumerate(fields):
            hit[f] = rec[f_i]

        if float(hit['i-Evalue']) <= dom_i_evalue_thresh:
            hits.append(hit)

    if hits:
        query_hits = hits
            
    return query_hits


# get_hmmer_dom_file_buf ()
#
def get_hmmer_dom_file_buf (file_path):
    outbuf = []

    with open (file_path, 'r') as f_h:
        for f_line in f_h:
            f_line = f_line.rstrip("\n")
            outbuf.append(f_line)
            
    return outbuf
    

# write_merged_tsv_outfile ()
#
def write_merged_tsv_outfile (outtsvfile, query_hits, fields):
    outbuf = []

    # header
    row = fields
    outbuf.append("\t".join(row))

    for hit in query_hits:
        #row = [gid]
        row = []
        for f in fields:
            row.append(hit[f])
        outbuf.append("\t".join(row))

                
    # write outfile
    with open(outtsvfile, 'w') as out_h:
        out_h.write("\n".join(outbuf)+"\n")
    
    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # config
    dom_i_evalue_thresh = .001
    fields = 'target_name target_accession tlen query_name query_accession qlen E-value full_bitscore full_bias dom_n dom_total c-Evalue i-Evalue dom_bitscore dom_bias hmm_from hmm_to ali_from ali_to env_from env_to acc description'.split()
    
    # read querys
    query_hits = read_query_hits (args.annotfile,
                                  dom_i_evalue_thresh,
                                  fields)

    # write merged tsv hit file
    out_file = write_merged_tsv_outfile (args.outtsvfile,
                                         query_hits,
                                         fields)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

