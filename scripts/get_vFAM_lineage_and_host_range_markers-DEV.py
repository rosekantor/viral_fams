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

    parser.add_argument("-f", "--faminfofile", help="viral fam info file: lineage and host range")
    parser.add_argument("-o", "--outtsvfile", help="output family lineage and host reange tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit (-1)

    if args.faminfofile is None:
        print ("must specify --{}\n".format('faminfofile'))
        args_pass = False
    elif not os.path.exists(args.faminfofile) or \
         not os.path.isfile(args.faminfofile) or \
         not os.path.getsize(args.faminfofile) > 0:
        print ("--{} {} must exist and not be empty\n".format('faminfofile', args.faminfofile))
        args_pass = False

    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_fam_info()
#
def read_fam_info (fam_info_file):
    fam_info = dict()
    fields = []
    
    with open (fam_info_file, 'r') as finfo_h:
        for finfo_line in finfo_h:
            finfo_line = finfo_line.strip()
            if not finfo_line or finfo_line.startswith('#'):
                continue
            if finfo_line.startswith("FAM\t"):
                fields = finfo_line.split("\t")
                continue
            rec = finfo_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_info_file))
            fam_id = rec[0]
            fam_info[fam_id] = dict()
            for fi,field in enumerate(fields):
                fam_info[fam_id][field] = rec[fi]
                
    return (fam_info, fields)


# get_lineage_fams()
#
def get_lineage_fams (fam_info, fields):
    lineage_fams = dict()

    for fam_id in sorted(fam_info.keys()):

        lca_node = fam_info[fam_id]['HITS_LCA_NODE']
        if lca_node not in lineage_fams:
            lineage_fams[lca_node] = [fam_id]
        else:
            lineage_fams[lca_node].append(fam_id)

    return lineage_fams


# write_lineage_fams_tsv_outfile ()
#
def write_lineage_fams_tsv_outfile (outtsvfile,
                                    lineage_fams):
    outbuf = []

    # get counts for sorting
    lineage_fam_counts = dict()
    for lca_node in lineage_fams.keys():
        lineage_fam_counts[lca_node] = len(lineage_fams[lca_node])

    # base_fields
    base_fields = ['FAM_COUNT', 'LINEAGE_NODE', 'FAMS']

#    # get host range cats
#    all_host_range_cats = dict()
#    for fam_id in fam_host_range_tally:
#        for host_range_cat in fam_host_range_tally[fam_id].keys():
#            all_host_range_cats[host_range_cat] = True
#    host_range_cats = sorted (all_host_range_cats.keys())

    all_fields = []
    all_fields.extend(base_fields)
#    all_fields.extend(host_range_cats)
    
    # header
    outbuf.append("\t".join(all_fields))

    # build output
    for lca_node in sorted(lineage_fam_counts, key=lineage_fam_counts.get, reverse=True):
        row = [str(lineage_fam_counts[lca_node])]
        row.append(lca_node)
        row.append(",".join(lineage_fams[lca_node]))
        outbuf.append("\t".join(row))
                
    # write outfile
    with open(outtsvfile, 'w') as out_h:
        out_h.write("\n".join(outbuf)+"\n")
    
    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # read fam info for lineage and host range and get stats
    (fam_info, fields) = read_fam_info (args.faminfofile)

    # get lineage fams
    lineage_fams = get_lineage_fams (fam_info, fields)

    # get host range fam tally
#    host_range_fam_tally = calc_host_range_tally (fam_info)
    
    # write lineage fams tsv outfile
    lineage_out_file = write_lineage_fams_tsv_outfile (args.outtsvfile,
                                                       lineage_fams)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

