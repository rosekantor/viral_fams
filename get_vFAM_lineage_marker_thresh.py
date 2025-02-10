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
    parser = argparse.ArgumentParser(description="determine bitscore thresholds for lineage markers")

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


# read_fam_genome_hits_info ()
#
def read_fam_genome_hits_info (fam_hits_file, coverage_thresh):
    fam_genome_hits_info = dict()
    fields = []
    #fi_lookup = dict()

    with open (fam_hits_file, 'r') as hmmer_h:
        for hmmer_line in hmmer_h:
            hmmer_line = hmmer_line.strip()
            if not hmmer_line or hmmer_line.startswith('#'):
                continue
            if hmmer_line.startswith('genome_id'):
                fields = hmmer_line.split("\t")
                #for fi,field in enumerate(fields):
                #    fi_lookup[field] = fi
                continue
            rec = hmmer_line.split("\t")
            if not fields:
                raise ValueError ("need fields header in file {}".format(fam_hits_file))
            this_hit = dict()
            for fi,field in enumerate(fields):
                this_hit[field] = rec[fi]

            fam_id = this_hit['target_name']
            genome_id = this_hit['genome_id']
            gene_id = this_hit['query_name']
            target_hmm_len = this_hit['tlen']
            target_hmm_ali_from = this_hit['hmm_from']
            target_hmm_ali_to = this_hit['hmm_to']

            if not valid_hmm_hit (coverage_thresh, target_hmm_len, target_hmm_ali_from, target_hmm_ali_to):
                continue
            
            if fam_id not in fam_genome_hits_info:
                fam_genome_hits_info[fam_id] = dict()
            if genome_id not in fam_genome_hits_info[fam_id]:
                fam_genome_hits_info[fam_id][genome_id] = []
                
            fam_genome_hits_info[fam_id][genome_id].append(this_hit)
            
    return (fam_genome_hits_info, fields)


# valid_hmm_hit()
#
def valid_hmm_hit (coverage_thresh, hmm_len, hmm_ali_from, hmm_ali_to):
    valid = False
    coverage_thresh = float(coverage_thresh)

    if float(int(hmm_ali_to)-int(hmm_ali_from)+1)/float(hmm_len) >= coverage_thresh:
        valid = True

    return valid

    
# calc_fam_bitscore_threshold_per_lineage_node ()
#
def calc_fam_bitscore_threshold_per_lineage_node (fam_genome_hits_info,
                                                  genome_info):
    fam_bs_vals_per_lineage_node = dict()
    fam_bs_thresh_per_lineage_node = dict()

    for fam_id in fam_genome_hits_info.keys():

        # get LCA base
        common_lineage = []
        for genome_id in fam_genome_hits_info[fam_id].keys():
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
            
        fam_lca = '; '.join(common_lineage)

        # store bitscore vals
        fam_bs_vals_per_lineage_node[fam_id] = dict()
        fam_bs_vals_per_lineage_node[fam_id][fam_lca] = dict()
        fam_bs_vals_per_lineage_node[fam_id][fam_lca]['vals'] = []
        child_start_i = len(common_lineage)

        for genome_id in fam_genome_hits_info[fam_id].keys():
            this_lineage = genome_info[genome_id]['v_lineage'].split('; ')


            # build structure
            for fam_hit in fam_genome_hits_info[fam_id][genome_id]:
                dom_bitscore = fam_hit['dom_bitscore']
                fam_bs_vals_per_lineage_node[fam_id][fam_lca]['vals'].append(dom_bitscore)

                parent_node = fam_bs_vals_per_lineage_node[fam_id][fam_lca]
                if child_start_i < len(this_lineage):
                    for taxon in this_lineage[child_start_i:]:
                        if 'child' not in parent_node:
                            parent_node['child'] = dict()
                        if taxon not in parent_node['child']:
                            parent_node['child'][taxon] = dict()
                            parent_node['child'][taxon]['vals'] = []
                        parent_node['child'][taxon]['vals'].append(dom_bitscore)
                        parent_node = parent_node['child'][taxon] 
                    
        # reduce scores to lowest per node and trim to main trunk
        this_fam_bs_per_lineage_node = min_val_bs_lineage (fam_bs_vals_per_lineage_node[fam_id][fam_lca])
        this_fam_bs_per_lineage_node = trim_bs_lineage (this_fam_bs_per_lineage_node)

        fam_bs_thresh_per_lineage_node[fam_id] = dict()
        fam_bs_thresh_per_lineage_node[fam_id][fam_lca] = this_fam_bs_per_lineage_node
        
    return fam_bs_thresh_per_lineage_node


# min_val_bs_lineage ()
#
def min_val_bs_lineage (fam_bs_vals_per_lineage_node):
    this_fam_bs_per_lineage_node = dict()

    this_fam_bs_per_lineage_node['val'] = min_float_strs(fam_bs_vals_per_lineage_node['vals'])

    if 'child' in fam_bs_vals_per_lineage_node:
        for taxon in fam_bs_vals_per_lineage_node['child'].keys():
            this_fam_bs_per_lineage_node['child'] = dict()
            this_fam_bs_per_lineage_node['child'][taxon] = min_val_bs_lineage (fam_bs_vals_per_lineage_node['child'][taxon]) 
        
    return this_fam_bs_per_lineage_node


# trim_bs_lineage ()
#
def trim_bs_lineage (fam_bs_per_lineage_node):
    this_trimmed_bs_lineage = dict()

    this_trimmed_bs_lineage['val'] = fam_bs_per_lineage_node['val']

    if 'child' in fam_bs_per_lineage_node:
        this_trimmed_bs_lineage['child'] = dict()
        children_taxa = fam_bs_per_lineage_node['child'].keys()
        best_child_taxa = None
        best_child_score = -100000000000000
        for taxon in children_taxa:
            this_score = fam_bs_per_lineage_node['child'][taxon]['val']
            if not best_child_taxa or this_score > best_child_score:
                best_child_taxa = [taxon]
                best_child_score = this_score
            elif this_score == best_child_score:
                best_child_taxa.append(taxon)
                
        for taxon in best_child_taxa:
            this_trimmed_bs_lineage['child'][taxon] = trim_bs_lineage(fam_bs_per_lineage_node['child'][taxon])

    return this_trimmed_bs_lineage
    

# min_float_strs ()
#
def min_float_strs (val_list):
    min_val = 1000000000000000000

    for val in val_list:
        if float(val) < min_val:
            min_val = float(val)
    return min_val        


# write_fam_bitscore_thresh_per_node ()
#
def write_fam_bitscore_thresh_per_node (outtsvfile, fam_bs_thresh_per_node):
    outbuf = []

    # header
    header = ['FAM', 'BITSCORE_THRESH', 'LINEAGE']
    outbuf.append("\t".join(header))

    # build output
    for fam_id in sorted(fam_bs_thresh_per_node.keys()):
        #print (fam_id, flush=True) # DEBUG
        taxon = list(fam_bs_thresh_per_node[fam_id].keys())[0]
        this_node = fam_bs_thresh_per_node[fam_id][taxon]

        bit_score = this_node['val']
        tax_string = taxon
        row = [fam_id, str(bit_score), tax_string]
        #print ("\t".join(row), flush=True)  # DEBUG
        outbuf.append("\t".join(row))

        last_bit_score = bit_score
        while 'child' in this_node:
            if len(this_node['child'].keys()) > 1:
                break
            taxon = list(this_node['child'].keys())[0]
            bit_score = this_node['child'][taxon]['val']
            if bit_score != last_bit_score:
                tax_string = tax_string+'; '+taxon
                row = [fam_id, str(bit_score), tax_string]
                outbuf.append("\t".join(row))
                #print ("\t".join(row), flush=True)  # DEBUG
            this_node = this_node['child'][taxon]
            last_bit_score = bit_score
            
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
    (fam_genome_hits_info, fam_hits_fields) = read_fam_genome_hits_info (args.famhitsfile,
                                                                         args.coverage)

    # get bitscore thresh for each fam at each lineage node
    fam_bs_thresh_per_node = calc_fam_bitscore_threshold_per_lineage_node (fam_genome_hits_info,
                                                                           genome_info)

    
    # write bitscore thresh for each fam per lineage
    out_file = write_fam_bitscore_thresh_per_node (args.outtsvfile,
                                                   fam_bs_thresh_per_node)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())

