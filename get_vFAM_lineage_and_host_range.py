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
    parser = argparse.ArgumentParser(description="determine lineage last-common-ancestor and host range hits for each viral protein family")

    parser.add_argument("-g", "--genomeinfofile", help="genome info file, lineage, host range, etc.")
    parser.add_argument("-f", "--famhitsfile", help="merged fam hmmer hits tsv table file")
    parser.add_argument("-c", "--coverage", help="coverage min (0.0-1.0, def=0.5)")
    parser.add_argument("-o", "--outtsvfile", help="output family lineage and host reange tsv table")
    
    args = parser.parse_args()
    args_pass = True

    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit (-1)

    if args.genomeinfofile is None:
        print ("must specify --{}\n".format('genomeinfofile'))
        args_pass = False
    elif not os.path.exists(args.genomeinfofile) or \
         not os.path.isfile(args.genomeinfofile) or \
         not os.path.getsize(args.genomeinfofile) > 0:
        print ("--{} {} must exist and not be empty\n".format('genomeinfofile', args.genomeinfofile))
        args_pass = False

    if args.famhitsfile is None:
        print ("must specify --{}\n".format('famhitsfile'))
        args_pass = False
    elif not os.path.exists(args.famhitsfile) or \
         not os.path.isfile(args.famhitsfile) or \
         not os.path.getsize(args.famhitsfile) > 0:
        print ("--{} {} must exist and not be empty\n".format('famhitsfile', args.famhitsfile))
        args_pass = False

    if args.outtsvfile is None:
        print ("must specify --{}\n".format('outtsvfile'))
        args_pass = False

    if args.coverage is None:
        args.coverage = '0.5'
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_genome_info()
#
def read_genome_info (genome_info_file):
    genome_info = dict()
    fields = []
    
    with open (genome_info_file, 'r') as ginfo_h:
        for ginfo_line in ginfo_h:
            ginfo_line = ginfo_line.strip()
            if not ginfo_line or ginfo_line.startswith('#'):
                continue
            if ginfo_line.startswith('rs_genome_id'):
                fields = ginfo_line.split("\t")
                continue
            rec = ginfo_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(genome_info_file))
            genome_id = rec[0]
            genome_info[genome_id] = dict()
            for fi,field in enumerate(fields):
                genome_info[genome_id][field] = rec[fi]
                
    return (genome_info, fields)


# read_fam_genome_hits ()
#
def read_fam_genome_hits (fam_hits_file, coverage_thresh):
    fam_genome_hits = dict()
    fields = []
    fi_lookup = dict()

    with open (fam_hits_file, 'r') as hmmer_h:
        for hmmer_line in hmmer_h:
            hmmer_line = hmmer_line.strip()
            if not hmmer_line or hmmer_line.startswith('#'):
                continue
            if hmmer_line.startswith('genome_id'):
                fields = hmmer_line.split("\t")
                for fi,field in enumerate(fields):
                    fi_lookup[field] = fi
                continue
            rec = hmmer_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_hits_file))
            fam_id = rec[fi_lookup['target_name']]
            genome_id = rec[fi_lookup['genome_id']]
            gene_id = rec[fi_lookup['query_name']]
            target_hmm_len = rec[fi_lookup['tlen']]
            target_hmm_ali_from = rec[fi_lookup['hmm_from']]
            target_hmm_ali_to = rec[fi_lookup['hmm_to']]

            if not valid_hmm_hit (coverage_thresh, target_hmm_len, target_hmm_ali_from, target_hmm_ali_to):
                continue
            
            if fam_id not in fam_genome_hits:
                fam_genome_hits[fam_id] = dict()
            if genome_id not in fam_genome_hits[fam_id]:
                fam_genome_hits[fam_id][genome_id] = dict()
                
            fam_genome_hits[fam_id][genome_id][gene_id] = True
            
    return fam_genome_hits


# valid_hmm_hit()
#
def valid_hmm_hit (coverage_thresh, hmm_len, hmm_ali_from, hmm_ali_to):
    valid = False
    coverage_thresh = float(coverage_thresh)

    if float(int(hmm_ali_to)-int(hmm_ali_from)+1)/float(hmm_len) >= coverage_thresh:
        valid = True

    return valid

    
# calc_lineage_last_common_ancestor_node()
#
def calc_lineage_last_common_ancestor_node (fam_genome_hits, genome_info):
    fam_lca = dict()

    for fam_id in fam_genome_hits.keys():
        common_lineage = []

        for genome_id in fam_genome_hits[fam_id].keys():
            this_lineage = genome_info[genome_id]['v_lineage'].split('; ')
            if not common_lineage:
                common_lineage.extend(this_lineage)
                continue
            new_common_lineage = []
            shorter_lineage_len = len(common_lineage)
            if len(this_lineage) < shorter_lineage_len:
                shorter_lineage_len = len(this_lineage) 
            for lineage_i in range(shorter_lineage_len):
                if this_lineage[lineage_i] == common_lineage[lineage_i]:
                    new_common_lineage.append(common_lineage[lineage_i])
                else:
                    break
            common_lineage = new_common_lineage
            
        fam_lca[fam_id] = '; '.join(common_lineage)

    return fam_lca


# calc_host_range_tally ()
#
def calc_host_range_tally (fam_genome_hits, genome_info):
    fam_host_range_tally = dict()

    for fam_id in fam_genome_hits.keys():
        fam_host_range_tally[fam_id] = dict()
        for genome_id in fam_genome_hits[fam_id].keys():
            this_host_range_cat = genome_info[genome_id]['man_h_category']

            if this_host_range_cat not in fam_host_range_tally[fam_id]:
                fam_host_range_tally[fam_id][this_host_range_cat] = 1
            else:
                fam_host_range_tally[fam_id][this_host_range_cat] += 1

    return fam_host_range_tally


# write_lineage_and_host_range_tsv_outfile ()
#
def write_lineage_and_host_range_tsv_outfile (outtsvfile,
                                              fam_genome_hits,
                                              fam_lca,
                                              fam_host_range_tally):
    outbuf = []

    # base_fields
    base_fields = ['FAM', 'HITS_LCA_NODE', 'GENOME_HITS']

    # get host range cats
    all_host_range_cats = dict()
    for fam_id in fam_host_range_tally:
        for host_range_cat in fam_host_range_tally[fam_id].keys():
            all_host_range_cats[host_range_cat] = True
    host_range_cats = sorted (all_host_range_cats.keys())

    all_fields = []
    all_fields.extend(base_fields)
    all_fields.extend(host_range_cats)
    
    # header
    outbuf.append("\t".join(all_fields))

    # build output
    for fam_id in sorted(fam_genome_hits.keys()):
        row = [fam_id]
        row.append(fam_lca[fam_id])

        genome_hits = []
        for genome_id in sorted(fam_genome_hits[fam_id].keys()):
            for gene_id in sorted(fam_genome_hits[fam_id][genome_id].keys()):
                genome_hits.append(":".join([genome_id,gene_id]))
        row.append(";".join(genome_hits))

        for host_range_cat in host_range_cats:
            if host_range_cat not in fam_host_range_tally[fam_id]:
                row.append('0')
            else:
                row.append(str(fam_host_range_tally[fam_id][host_range_cat]))

        outbuf.append("\t".join(row))
                
    # write outfile
    with open(outtsvfile, 'w') as out_h:
        out_h.write("\n".join(outbuf)+"\n")
    
    return outtsvfile


# main()
#
def main() -> int:
    args = getargs()

    # get genomeids
    (genome_info, genome_fields) = read_genome_info (args.genomeinfofile)

    # read fam hits to genomes and genes
    fam_genome_hits = read_fam_genome_hits (args.famhitsfile, args.coverage)

    # get lineage last common ancestor
    fam_lca = calc_lineage_last_common_ancestor_node (fam_genome_hits,
                                                      genome_info)

    # get host range cat tally
    fam_host_range_tally = calc_host_range_tally (fam_genome_hits,
                                                  genome_info)
    
    # write lineage and host range tsv outfile
    out_file = write_lineage_and_host_range_tsv_outfile (args.outtsvfile,
                                                         fam_genome_hits,
                                                         fam_lca,
                                                         fam_host_range_tally)
    
    #print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

