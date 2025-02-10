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

    parser.add_argument("-t", "--targetfam", help="viral target fam file")
    parser.add_argument("-f", "--famhitfile", help="viral fam hit file")
    parser.add_argument("-p", "--pdbhitfile", help="pdb hit file")
    parser.add_argument("-o", "--outtsvfile", help="output fam and pdb hit tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 9:
        parser.print_help()
        sys.exit (-1)

    if args.targetfam is None:
        print ("must specify --{}\n".format('targetfam'))
        args_pass = False

    if args.famhitfile is None:
        print ("must specify --{}\n".format('famhitfile'))
        args_pass = False
    elif not os.path.exists(args.famhitfile) or \
         not os.path.isfile(args.famhitfile) or \
         not os.path.getsize(args.famhitfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('famhitfile', args.famhitfile))
        args_pass = False

    if args.pdbhitfile is None:
        print ("must specify --{}\n".format('pdbhitfile'))
        args_pass = False
    elif not os.path.exists(args.pdbhitfile) or \
         not os.path.isfile(args.pdbhitfile) or \
         not os.path.getsize(args.famhitfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('pdbhitfile', args.pdbhitfile))
        args_pass = False

    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_target_fam_gene_hits ()
#
def read_target_fam_gene_hits (target_fam, fam_hit_file):
    fam_gene_hits = dict()
    fields = []
    
    with open (fam_hit_file, 'r') as fhit_h:
        for fhit_line in fhit_h:
            fhit_line = fhit_line.strip()
            if not fhit_line or fhit_line.startswith('#'):
                continue
            if fhit_line.startswith("genome_id\t"):
                fields = fhit_line.split("\t")
                continue
            rec = fhit_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_hit_file))
            this_hit = dict()
            for fi,field in enumerate(fields):
                this_hit[field] = rec[fi]
            genome_id = this_hit['genome_id']
            gene_id = this_hit['query_name']
            fam_id = this_hit['target_name']

            if fam_id != target_fam:
                continue

            if genome_id not in fam_gene_hits:
                fam_gene_hits[genome_id] = dict()
            if gene_id not in fam_gene_hits[genome_id]:
                fam_gene_hits[genome_id][gene_id] = []
            fam_gene_hits[genome_id][gene_id].append(this_hit)
            
    return (fam_gene_hits, fields)


# read_pdb_all_hits()
#
def read_pdb_all_hits (pdb_hit_file, fam_gene_hits, fam_hits_fields):
    pdb_all_hits = dict()
    fields = []
    
    with open (pdb_hit_file, 'r') as phit_h:
        for phit_line in phit_h:
            phit_line = phit_line.strip()
            if not phit_line or phit_line.startswith('#'):
                continue
            if phit_line.startswith("v_genome_id\t"):
                fields = phit_line.split("\t")
                continue
            rec = phit_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(pdb_hit_file))
            this_hit = dict()
            for fi,field in enumerate(fields):
                this_hit[field] = rec[fi]
            v_genome_id = this_hit['v_genome_id']
            v_gene_id = this_hit['v_gene_id']
            pdb_hit_id = this_hit['pdb_hit_id']

            if v_genome_id not in fam_gene_hits:
                continue
            if v_gene_id not in fam_gene_hits[v_genome_id]:
                continue

            for fam_hit in fam_gene_hits[v_genome_id][v_gene_id]:
                if check_overlap (0.5,
                                  fam_hit['env_from'],
                                  fam_hit['env_to'],
                                  this_hit['pdb_hit_q_beg'],
                                  this_hit['pdb_hit_q_end']):

                    pdb_base = pdb_hit_id[:4]
                    
                    if v_genome_id not in pdb_all_hits:
                        pdb_all_hits[v_genome_id] = dict()
                    if v_gene_id not in pdb_all_hits[v_genome_id]:
                        pdb_all_hits[v_genome_id][v_gene_id] = []

                    # fix chain ids    
                    this_hit['pdb_hit_id'] = fix_pdb_id (this_hit['pdb_hit_id'])
                    this_hit['cochain_pdb_id'] = fix_pdb_id (this_hit['cochain_pdb_id'])

                    # store hit
                    new_hit = dict()
                    for field in fields:
                        new_hit[field] = this_hit[field]
                    
                    for fam_hit_field in fam_hits_fields:
                        if fam_hit_field in new_hit:
                            raise ValueError ("overlapping field name {}".format(fam_hit_field))
                        new_hit[fam_hit_field] = fam_hit[fam_hit_field]
                    pdb_all_hits[v_genome_id][v_gene_id].append(new_hit)
            
    return (pdb_all_hits, fields)


# check_overlap ()
#
def check_overlap (overlap_min, region_1_beg, region_1_end, region_2_beg, region_2_end):
    overlap_met = False
 
    region_1_beg = int(region_1_beg)
    region_1_end = int(region_1_end)
    region_2_beg = int(region_2_beg)
    region_2_end = int(region_2_end)
   
    region_1_len = region_1_end - region_1_beg + 1
    region_2_len = region_2_end - region_2_beg + 1

    shorter_len = region_1_len
    if region_2_len < region_1_len:
        shorter_len = region_2_len

    if region_2_end < region_1_beg:
        return False
    if region_2_beg > region_1_end:
        return False

    if region_2_end >= region_1_beg and region_2_end <= region_1_end:
        if region_2_beg >= region_1_beg:
            overlap_len = region_2_len
        else:
            overlap_len = region_2_end - region_1_beg + 1
    elif region_1_end >= region_2_beg:
        if region_1_beg >= region_2_beg:
            overlap_len = region_1_len
        else:
            overlap_len = region_1_end - region_2_beg + 1

    if float(overlap_len) / float(shorter_len) >= overlap_min:
        overlap_met = True

    return overlap_met


# fix_pdb_id ()
#
def fix_pdb_id (pdb_id):
    new_pdb_id = None

    split_pdb_id = pdb_id.split('_')
    new_pdb_id = split_pdb_id[0]+'_'+','.join(list(split_pdb_id[1]))

    return new_pdb_id


# write_fam_hit_and_pdb_tsv_outfile ()
#
def write_fam_hit_and_pdb_tsv_outfile (outtsvfile,
                                       pdb_all_hits,
                                       fam_gene_hits_fields,
                                       pdb_hits_fields):
    outbuf = []

    # header
    header_row = fam_gene_hits_fields + pdb_hits_fields
    outbuf.append("\t".join(header_row))
    
    # genes
    for genome_id in sorted(pdb_all_hits.keys()):
        for gene_id in sorted(pdb_all_hits[genome_id].keys()):
            for hit in pdb_all_hits[genome_id][gene_id]:

                out_row = []
                for fam_gene_hit_field in fam_gene_hits_fields:
                    out_row.append(hit[fam_gene_hit_field])
                for pdb_hit_field in pdb_hits_fields:
                    out_row.append(hit[pdb_hit_field])
                    
                outbuf.append("\t".join(out_row))


    # write outfile
    if not outtsvfile:
        outtsvfile = 'STDOUT'
        print ("\n".join(outbuf))
    else:
        with open(outtsvfile, 'w') as out_h:
            out_h.write("\n".join(outbuf)+"\n")

    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # fam gene hits
    (fam_gene_hits, fam_gene_hits_fields) = read_target_fam_gene_hits (args.targetfam,
                                                                       args.famhitfile)

    # read top hit to pdb and store with fam_gene_hits info
    (pdb_all_hits, pdb_hits_fields) = read_pdb_all_hits (args.pdbhitfile,
                                                         fam_gene_hits,
                                                         fam_gene_hits_fields)

    # write fam and pdb hits tsv outfile
    out_file = write_fam_hit_and_pdb_tsv_outfile (args.outtsvfile,
                                                  pdb_all_hits,
                                                  fam_gene_hits_fields,
                                                  pdb_hits_fields)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

